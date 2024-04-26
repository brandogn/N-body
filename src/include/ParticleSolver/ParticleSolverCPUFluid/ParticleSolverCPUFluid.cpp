#include "ParticleSolverCPUFluid.h"
#include <glm/gtx/norm.hpp>
#include <cmath>

ParticleSolverCPUFluid::ParticleSolverCPUFluid(GridCPU *grid, float stepSize,
                                               float squaredSoft)
    : ParticleSolver() {
  this->squaredSoftening = squaredSoft;
  this->timeStep = stepSize;
  this->G = 1.0f;
  this->smoothingRadius = 0.2f;
  this->grid = grid;
}

void ParticleSolverCPUFluid::updateParticlePositions(
    ParticleSystem *particles) {
  this->grid->updateGrid(particles);

  #pragma omp parallel for schedule(static) shared(particles)
  for(size_t i =  0; i < particles->size(); i++){
    this->computeDensityMap(particles, i);
    this->computeGravityForce(particles, i);
  }

  #pragma omp parallel for schedule(static) shared(particles)
  for(size_t i =  0; i<particles->size(); i++){
    particles->updateParticlePosition(i, this->timeStep);
  }
}

float ParticleSolverCPUFluid::densityKernel(float distance, float radius) {
  if (distance < radius) {
    float scale = 15 / (2 * PI * std::pow(std::abs(radius), 5));
    float v = radius - distance;
    return v * v * scale;
  }
  return 0;
}

float ParticleSolverCPUFluid::densityDerivative(float distance, float radius) {
  if (distance <= radius) {
    float scale = 15 / (PI * std::pow(std::abs(radius), 5));
    float v = radius - distance;
    return -v * scale;
  }
  return 0;
}

float ParticleSolverCPUFluid::nearDensityKernel(float distance, float radius) {
  if (distance < radius) {
    float scale = 15 / (PI * std::pow(std::abs(radius), 6));
    float v = radius - distance;
    return v * v * v * scale;
  }
  return 0;
}

float ParticleSolverCPUFluid::nearDensityDerivative(float distance, float radius) {
  if (distance <= radius) {
    float scale = 45 / (PI * std::pow(std::abs(radius), 6));
    float v = radius - distance;
    return -v * v * scale;
  }
  return 0;
}

void ParticleSolverCPUFluid::computeDensityMap(ParticleSystem *particles, const unsigned int particleID) {
  glm::vec4 particlePosition = particles->getPositions()[particleID];
  
  // Get the bucket where the particle is located
  Bucket* bucket = this->grid->getBucketByPosition(particlePosition);

  // Compute the density contribution from neighboring particles in the bucket for each particle 
  float density = 0.f;
  float near_density = 0.f;
  for (size_t j = 0; j < bucket->getNumParticles(); j++) {
    const unsigned int otherParticleId = bucket->getParticleId(j);
    const glm::vec4 vector_i_j = particles->getPositions()[otherParticleId] - particlePosition;
    const float distance_i_j = std::pow(glm::length2(vector_i_j) + this->squaredSoftening, 1.5);
    density += densityKernel(distance_i_j, smoothingRadius);
    near_density += nearDensityKernel(distance_i_j, smoothingRadius);
  }

  particles->getDensities()[particleID] = glm::vec4(density, near_density, 0.f, 0.f);
}



void ParticleSolverCPUFluid::computeGravityForce(
    ParticleSystem *particles, const unsigned int particleId) {

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

  particles->getForces()[particleId] = totalForce;
}

bool ParticleSolverCPUFluid::usesGPU() { return false; }

ParticleSolverCPUFluid::~ParticleSolverCPUFluid() noexcept {
    delete this->grid;
}

float ParticleSolverCPUFluid::getSquaredSoftening() {
  return this->squaredSoftening;
}