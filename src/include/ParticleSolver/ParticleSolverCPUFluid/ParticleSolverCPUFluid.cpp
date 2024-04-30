#include "ParticleSolverCPUFluid.h"
#include <glm/gtx/norm.hpp>
#include <cmath>
#include <iostream>

ParticleSolverCPUFluid::ParticleSolverCPUFluid(GridCPU *grid, float stepSize,
                                               float squaredSoft)
    : ParticleSolver() {
  this->squaredSoftening = squaredSoft;
  this->timeStep = stepSize;
  this->G = 1.0f;
  this->smoothingRadius = 0.1f;
  this->targetDensity = 20.f;
  this->pressureMultiplier = 30.f;
  this->nearPressureMultiplier = 1.95f;
  this->viscosityStrength = 0.f;
  this->grid = grid;
}

void ParticleSolverCPUFluid::updateParticlePositions(
    ParticleSystem *particles) {
  this->grid->updateGrid(particles);

  #pragma omp parallel for schedule(static) shared(particles)
  for(size_t i =  0; i < particles->size(); i++){
    this->computeTemperatures(particles, i);
    this->computeDensityMap(particles, i);
  }

  #pragma omp parallel for schedule(static) shared(particles)
  for(size_t i =  0; i < particles->size(); i++){
    this->computePressureForce(particles, i);
    this->computeViscosityForce(particles, i);
    this->computeGravityForce(particles, i);
    // std::cout << "DEBUG: {" 
    //   << particles->getForces()[i].x 
    //   << particles->getForces()[i].y 
    //   << particles->getForces()[i].z 
    //   << "} \n";
  }

  #pragma omp parallel for schedule(static) shared(particles)
  for(size_t i =  0; i<particles->size(); i++){
    particles->updateParticlePosition(i, this->timeStep);
  }
}

// void ParticleSolverCPUFluid::initializeShockWave(ParticleSystem *particles) {
//     float shockRadius = 0.5f; // Radius from center to differentiate between the inner and outer layers
//     float shockStrength = 100.0f; // Initial outward velocity for outer particles

//     for (size_t i = 0; i < particles->size(); i++) {
//         glm::vec4 particlePosition = particles->getPositions()[i];
//         float distanceToCenter = glm::length(particlePosition - glm::vec4(0.f));
//         if (distanceToCenter > shockRadius) {
//             glm::vec4 direction = glm::normalize(particlePosition - glm::vec4(0.f));
//             particles->getVelocities()[i] += direction * shockStrength; // Impart initial outward velocity
//         }
//     }
// }

float ParticleSolverCPUFluid::pressureFromDensity(float density) {
  return (density - targetDensity) * pressureMultiplier;
}

float ParticleSolverCPUFluid::nearPressureFromDensity(float nearDensity) {
  return nearDensity * nearPressureMultiplier;
}

void ParticleSolverCPUFluid::computeDensityMap(ParticleSystem *particles, const unsigned int particleID) {
  glm::vec4 particlePosition = particles->getPositions()[particleID];
  
  // Get the bucket where the particle is located
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);

  // Compute the density contribution from neighboring particles 
  float density = 0.f;
  float near_density = 0.f;
  for (size_t j = 0; j < bucket->getNumParticles(); j++) {
    const unsigned int otherParticleId = bucket->getParticleId(j);
    const float distance_i_j = glm::distance(particles->getPositions()[otherParticleId], particlePosition);
    if (distance_i_j <= smoothingRadius) {
      density += FluidMath::densityKernel(distance_i_j, smoothingRadius);
      near_density += FluidMath::nearDensityKernel(distance_i_j, smoothingRadius);
    }
  }

  particles->getDensities()[particleID] = glm::vec4(density, near_density, 0.f, 0.f);
}

void ParticleSolverCPUFluid::computeTemperatures(ParticleSystem *particles, const unsigned int particleId) {
  glm::vec4 particlePosition = particles->getPositions()[particleId];

  // Get the bucket where the particle is located
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);

  for (size_t j = 0; j < bucket->getNumParticles(); j++) {
    particles->getTemperatures()[particleId].x = FluidMath::computeTemperature(particles->getMasses()[particleId].x, particles->getVelocities()[particleId].x);
  }
}

void ParticleSolverCPUFluid::computePressureForce(ParticleSystem *particles, const unsigned int particleID) {

  /*
  glm::vec4 particlePosition = particles->getPositions()[particleID];
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);


  float distanceToCenter = glm::length(particlePosition - glm::vec4(0.f)); // is center of the sim at origin?
  const float thresholdRadius = 0.5f; // adjust this after running, whatever looks more right - 0.5 is default
  float explosionPressureMultiplier = (distanceToCenter > thresholdRadius) ? 10.0f : 1.0f;

  const float density = particles->getDensities()[particleID].x;
  const float pressure = explosionPressureMultiplier * pressureFromDensity(density);
  glm::vec4 pressureForce(0.f);

  
  for (size_t j = 0; j < bucket->getNumParticles(); j++) {
    const unsigned int otherParticleId = bucket->getParticleId(j);
    if (otherParticleId != particleID) {
      glm::vec4 vector_i_j = particles->getPositions()[otherParticleId] - particlePosition;
      const float distance_i_j = glm::distance(particles->getPositions()[otherParticleId], particlePosition);

      if (distance_i_j < smoothingRadius) {
        const float densityNeighbor = particles->getDensities()[otherParticleId].x;
        const float nearDensityNeighbor = particles->getDensities()[otherParticleId].y;
        const float neighborPressure = pressureFromDensity(densityNeighbor);
        const float neighborNearPressure = nearPressureFromDensity(nearDensityNeighbor);

        const float sharedPressure = (pressure + neighborPressure) / 2;
        const float sharedNearPressure = (nearPressure + neighborNearPressure) / 2;

        if (distance_i_j <= 0) {
          vector_i_j = glm::vec4(0.f, 1.f, 0.f, 0.f);
        }

        pressureForce += vector_i_j * FluidMath::densityDerivative(distance_i_j, smoothingRadius) * sharedPressure / densityNeighbor;
        pressureForce += vector_i_j * FluidMath::nearDensityDerivative(distance_i_j, smoothingRadius) * sharedNearPressure / nearDensityNeighbor;
      }
    }
  }

  glm::vec4 acceleration = pressureForce / density;
  particles->getVelocities()[particleID] += acceleration * this->timeStep;
  */


  glm::vec4 particlePosition = particles->getPositions()[particleID];
  
  // Get the bucket where the particle is located
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);

  // Compute the pressure forces exerted by neighboring particles
  const float density = particles->getDensities()[particleID].x;
  const float nearDensity = particles->getDensities()[particleID].y;
  const float pressure = pressureFromDensity(density);
  const float nearPressure = nearPressureFromDensity(nearDensity);
  glm::vec4 pressureForce (0.f);

  for (size_t j = 0; j < bucket->getNumParticles(); j++) {
    const unsigned int otherParticleId = bucket->getParticleId(j);
    if (otherParticleId != particleID) {
      glm::vec4 vector_i_j = particles->getPositions()[otherParticleId] - particlePosition;
      const float distance_i_j = glm::distance(particles->getPositions()[otherParticleId], particlePosition);

      if (distance_i_j <= smoothingRadius) {
        const float densityNeighbor = particles->getDensities()[otherParticleId].x;
        const float nearDensityNeighbor = particles->getDensities()[otherParticleId].y;
        const float neighborPressure = pressureFromDensity(densityNeighbor);
        const float neighborNearPressure = nearPressureFromDensity(nearDensityNeighbor);

        const float sharedPressure = (pressure + neighborPressure) / 2;
        const float sharedNearPressure = (nearPressure + neighborNearPressure) / 2;

        if (distance_i_j <= 0) {
          vector_i_j = glm::vec4(0.f, 1.f, 0.f, 0.f);
        }

        pressureForce += vector_i_j * FluidMath::densityDerivative(distance_i_j, smoothingRadius) * sharedPressure / densityNeighbor;
        pressureForce += vector_i_j * FluidMath::nearDensityDerivative(distance_i_j, smoothingRadius) * sharedNearPressure / nearDensityNeighbor;
      }
    }
  }

  glm::vec4 acceleration = pressureForce / density;
  particles->getVelocities()[particleID] += acceleration * this->timeStep;
}

void ParticleSolverCPUFluid::computeViscosityForce(ParticleSystem *particles, const unsigned int particleId) {
  glm::vec4 particlePosition = particles->getPositions()[particleId];
  glm::vec4 viscosityForce (0.f);

  glm::vec4 particleVelocity = particles->getVelocities()[particleId];
  
  // Get the bucket where the particle is located
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);

  for (size_t j = 0; j < bucket->getNumParticles(); j++) {
    const unsigned int otherParticleId = bucket->getParticleId(j);
    if (otherParticleId != particleId) {
      glm::vec4 vector_i_j = particles->getPositions()[otherParticleId] - particlePosition;
      const float distance_i_j = glm::distance(particles->getPositions()[otherParticleId], particlePosition);

      if (distance_i_j < smoothingRadius) {
        const glm::vec4 neighborVelocity = particles->getVelocities()[otherParticleId];
        viscosityForce += (neighborVelocity - particleVelocity) * FluidMath::viscosityKernel(distance_i_j, smoothingRadius);
      }
    }
  }
  particles->getVelocities()[particleId] += viscosityForce * viscosityStrength * this->timeStep;
}

void ParticleSolverCPUFluid::computeGravityForce(
    ParticleSystem *particles, const unsigned int particleId) {


  /*
  glm::vec4 particlePosition = particles->getPositions()[particleId];
  glm::vec4 centerOfMass(0.f, 0.f, 0.f, 0.f); // again assuming center of mass at origin

  float distanceToCenter = glm::length(particlePosition - centerOfMass);
  float collapseThresholdRadius = 0.5f; // adjust if necessary

  // gravitational constants for normal and 'enhanced' gravitational pull near center
  float normalGravitationalStrength = G;
  float enhancedGravitationalStrength = G * 10.0f; // stronger gravitational pull within the collapse threshold, again can adjust this

  float gravitationalStrength = (distanceToCenter <= collapseThresholdRadius) ? enhancedGravitationalStrength : normalGravitationalStrength;

  glm::vec4 gravityDirection = glm::normalize(centerOfMass - particlePosition);
  glm::vec4 gravityForce = gravityDirection * gravitationalStrength * particles->getMasses()[particleId];

  // apply central gravity force towards the simulation center
  particles->getForces()[particleId] += gravityForce;

  // tried to follow same format as u guys for rest of this
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);
    
  // inside bucket
  for(size_t j = 0; j < bucket->getNumParticles(); j++){
      const unsigned int otherParticleId = bucket->getParticleId(j);
      if (otherParticleId != particleId) {
          glm::vec4 otherParticlePosition = particles->getPositions()[otherParticleId];
          glm::vec4 vector_i_j = otherParticlePosition - particlePosition;
          float distance_i_j = glm::length(vector_i_j);

          if (distance_i_j > 0) {
              float forceMagnitude = (normalGravitationalStrength * particles->getMasses()[particleId] * particles->getMasses()[otherParticleId]) /
                                       (distance_i_j * distance_i_j + this->squaredSoftening);
              glm::vec4 forceVector = glm::normalize(vector_i_j) * forceMagnitude;
              particles->getForces()[particleId] -= forceVector;
              particles->getForces()[otherParticleId] += forceVector;
          }
      }
  }

  // other buckets
  for(size_t bucketId = 0; bucketId < this->grid->getTotalBuckets(); bucketId++){
      Bucket *otherBucket = this->grid->getBucketById(bucketId);
      if (bucket->getBucketId() != otherBucket->getBucketId()){
          for(size_t k = 0; k < otherBucket->getNumParticles(); k++) {
              unsigned int otherParticleId = otherBucket->getParticleId(k);
              glm::vec4 otherParticlePosition = particles->getPositions()[otherParticleId];
              glm::vec4 vector_i_j = otherParticlePosition - particlePosition;
              float distance_i_j = glm::length(vector_i_j);

              if (distance_i_j > 0) {
                  float forceMagnitude = (normalGravitationalStrength * particles->getMasses()[particleId] * particles->getMasses()[otherParticleId]) /
                                           (distance_i_j * distance_i_j + this->squaredSoftening);
                  glm::vec4 forceVector = glm::normalize(vector_i_j) * forceMagnitude;
                  particles->getForces()[particleId] -= forceVector;
                  particles->getForces()[otherParticleId] += forceVector;
              }
          }
      }
  }
  */

  glm::vec4 particlePosition = particles->getPositions()[particleId];
  glm::vec4 totalForce (0.f);

  // Get the bucket where the particle is located
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);

  // Compute forces inside the bucket
  for(size_t j = 0; j < bucket->getNumParticles(); j++){
    const unsigned int otherParticleId = bucket->getParticleId(j);
    const glm::vec4 vector_i_j = particles->getPositions()[otherParticleId] - particlePosition;
    const float distance_i_j = std::pow(glm::length2(vector_i_j) + this->squaredSoftening, 1.5);
    totalForce += ((G * particles->getMasses()[otherParticleId].x) / distance_i_j) * vector_i_j;
  }

  // Compute the forces with other buckets
  Bucket *otherBucket = nullptr;
  for(size_t bucketId = 0; bucketId < this->grid->getTotalBuckets(); bucketId++){
    otherBucket = this->grid->getBucketById(bucketId);
    if (bucket->getBucketId() != otherBucket->getBucketId()){
      const glm::vec4 centerOfMass = otherBucket->getCenterOfMass();
      const float mass = centerOfMass.w;
      const glm::vec4 centerOfMassPosition = glm::vec4(centerOfMass.x, centerOfMass.y, centerOfMass.z, 0.f);
      const glm::vec4 vector_i_j = centerOfMassPosition - particlePosition;
      const float distance_i_j = std::pow(glm::length2(vector_i_j) + this->squaredSoftening, 1.5);
      totalForce += ((G * mass) / distance_i_j) * vector_i_j;
    }
  }

  particles->getForces()[particleId] += totalForce;
}

bool ParticleSolverCPUFluid::usesGPU() { return false; }

ParticleSolverCPUFluid::~ParticleSolverCPUFluid() noexcept {
    delete this->grid;
}

float ParticleSolverCPUFluid::getSquaredSoftening() {
  return this->squaredSoftening;
}
