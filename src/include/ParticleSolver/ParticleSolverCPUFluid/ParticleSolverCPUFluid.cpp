
#include "ParticleSolverCPUFluid.h"
#include <glm/gtx/norm.hpp>

ParticleSolverCPUFluid::ParticleSolverCPUFluid(float stepSize,
                                               float squaredSoft)
    : ParticleSolver() {
  this->squaredSoftening = squaredSoft;
  this->timeStep = stepSize;
  this->G = 1.0f;
}

void ParticleSolverCPUFluid::updateParticlePositions(
    ParticleSystem *particles) {
  for (size_t i = 0; i < particles->size(); i++) {
    this->computeGravityForce(particles, i);
  }
  for (size_t i = 0; i < particles->size(); i++) {
    particles->updateParticlePosition(i, this->timeStep);
  }
}

void ParticleSolverCPUFluid::computeDensityMap(ParticleSystem *particles) {}

void ParticleSolverCPUFluid::computeGravityForce(
    ParticleSystem *particles, const unsigned int particleId) {

  glm::vec4 particlePosition = particles->getPositions()[particleId];

  glm::vec4 totalForce(0.f);

  // for(size_t j = 0; j < particles->size(); j++){
  //     const glm::vec4 vector_i_j = particles->getPositions()[j] -
  //     particlePosition; const float distance_i_j =
  //     std::pow(glm::length2(vector_i_j) + this->squaredSoftening, 1.5);
  //     totalForce += ((G * particles->getMasses()[j].x) / distance_i_j) *
  //     vector_i_j;
  // }

  totalForce = particles->getVelocities()[particleId] * -1.0f;
  particles->getForces()[particleId] = totalForce;
}

bool ParticleSolverCPUFluid::usesGPU() { return false; }

float ParticleSolverCPUFluid::getSquaredSoftening() {
  return this->squaredSoftening;
}