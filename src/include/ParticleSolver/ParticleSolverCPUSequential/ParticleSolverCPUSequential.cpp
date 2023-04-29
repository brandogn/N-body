
#include "ParticleSolverCPUSequential.h"
#include <glm/gtx/norm.hpp>

ParticleSolverCPUSequential::ParticleSolverCPUSequential(float stepSize, float squaredSoft): ParticleSolver() {
    this->squaredSoftening = squaredSoft;
    this->timeStep = stepSize;
    this->G = 1.0f;
}

void ParticleSolverCPUSequential::updateParticlePositions(ParticleSystem *particles){
    for(size_t i = 0; i<particles->size(); i++){
        particles->updateParticlePosition(i, this->timeStep, this->computeGravityForce(particles, i));
    }
}

glm::vec4
ParticleSolverCPUSequential::computeGravityForce(ParticleSystem *particles, const unsigned int particleId) {

    glm::vec4 particlePosition = particles->getPositions()[particleId];
    float particleMass = particles->getMasses()[particleId].x;

    glm::vec4 totalForce (0.f);

    for(size_t j = 0; j < particles->size(); j++){
        const glm::vec4 vector_i_j = particles->getPositions()[j] - particlePosition;
        const float distance_i_j = std::pow(glm::length2(vector_i_j) + this->squaredSoftening, 1.5);
        totalForce += ((G * particleMass * particles->getMasses()[j].x) / distance_i_j) * vector_i_j;
    }

    return totalForce;
}

bool ParticleSolverCPUSequential::usesGPU() {return false;}


float ParticleSolverCPUSequential::getSquaredSoftening() {
    return this->squaredSoftening;
}