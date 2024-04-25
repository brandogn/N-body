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
    this->computeGravityForce(particles, i);
  }

  #pragma omp parallel for schedule(static) shared(particles)
  for(size_t i =  0; i<particles->size(); i++){
    particles->updateParticlePosition(i, this->timeStep);
  }
}

float ParticleSolverCPUFluid::smoothingKernelDensity(float distance, float radius) {
  if (distance < radius) {
    float scale = 315 / (64 * 3.14159265358979323846 * std::pow(std::abs(radius), 9));
    float v = radius * radius - distance - distance;
    return v * v * v * scale;
  }
  return 0;
}

void ParticleSolverCPUFluid::computeDensityMap(ParticleSystem *particles) {

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