
#include "ParticleSolverGPU.h"
#include <glad/glad.h>

ParticleSolverGPU::ParticleSolverGPU(float stepSize, float squaredSoft, std::string &pathToComputeShader): ParticleSolver() {
    this->computeShader = new ComputeShader(pathToComputeShader);
    this->computeShader->use();
    this->computeShader->setFloat("deltaTime", stepSize);
    this->computeShader->setFloat("squaredSoftening", squaredSoft);
}

void ParticleSolverGPU::updateParticlePositions(ParticleSystem *particles) {
    this->computeShader->use();
    this->computeShader->setInt("numParticles", particles->size());
    glDispatchCompute(ceil(particles->size() / 64.0), 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}


ParticleSolverGPU::~ParticleSolverGPU() noexcept {
    delete this->computeShader;
}

bool ParticleSolverGPU::usesGPU() {return true;}

float ParticleSolverGPU::getSquaredSoftening() {
    return this->squaredSoftening;
}