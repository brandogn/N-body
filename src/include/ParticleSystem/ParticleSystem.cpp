#include "ParticleSystem.h"

ParticleSystem::ParticleSystem(std::vector<Particle> &particles) {
    this->numParticles = particles.size();
    this->velocities = new glm::vec4[this->numParticles]();
    this->accelerations = new glm::vec4[this->numParticles]();
    this->positions = new glm::vec4[this->numParticles]();
    this->masses = new glm::vec4[this->numParticles]();
    this->densities = new glm::vec4[this->numParticles]();
    this->temperatures = new glm::vec4[this->numParticles]();
    this->forces = new glm::vec4[this->numParticles]();
    
    for (int i = 0; i < this->numParticles; i++) {
        this->velocities[i] = particles[i].velocity;
        this->accelerations[i] = particles[i].acceleration;
        this->positions[i] = particles[i].position;
        this->masses[i] = glm::vec4(particles[i].mass, 0, 0, 0);
        this->densities[i] = glm::vec4(0.f);
        this->temperatures[i] = glm::vec4(0.f);
        this->forces[i] = glm::vec4(0.f);
    }
}

ParticleSystem::ParticleSystem(ParticleSystem * other) {
    this->numParticles = other->size();
    this->velocities = new glm::vec4[this->numParticles]();
    this->accelerations = new glm::vec4[this->numParticles]();
    this->positions = new glm::vec4[this->numParticles]();
    this->masses = new glm::vec4[this->numParticles]();
    this->densities = new glm::vec4[this->numParticles]();
    this->temperatures = new glm::vec4[this->numParticles]();
    this->forces = new glm::vec4[this->numParticles]();

    for (int i = 0; i < this->numParticles; i++) {
        this->velocities[i] = other->getVelocities()[i];
        this->accelerations[i] = other->getAccelerations()[i];
        this->positions[i] = other->positions[i];
        this->masses[i] = other->masses[i];
        this->densities[i] = other->densities[i];
        this->temperatures[i] = other->temperatures[i];
        this->forces[i] = other->forces[i];
    }
}



/**
 * Updates a particle position
 * Performs the leapfrog integration
 * @param particleId
 * @param deltaTime
 */
void ParticleSystem::updateParticlePosition(unsigned int particleId, float deltaTime) {
    float dtDividedBy2 = deltaTime * 0.5f;

    // Compute velocity (i + 1/2)
    this->velocities[particleId] += this->accelerations[particleId] * dtDividedBy2;

    // Compute next position (i+1)
    this->positions[particleId] += this->velocities[particleId] * deltaTime;

    // Update acceleration (i+1)
    // F = MA;
    // A = F/M;  M is cancelled when calculating gravity force
    this->accelerations[particleId] = this->forces[particleId];

    // Compute next velocity (i+1)
    this->velocities[particleId] += this->accelerations[particleId] * dtDividedBy2;

}


unsigned int ParticleSystem::size() const {
    return this->numParticles;
}

glm::vec4* ParticleSystem::getMasses(){
    return this->masses;
}

glm::vec4* ParticleSystem::getPositions() {
    return this->positions;
}

glm::vec4* ParticleSystem::getVelocities() {
    return this->velocities;
}

glm::vec4* ParticleSystem::getAccelerations() {
    return this->accelerations;
}

glm::vec4* ParticleSystem::getDensities() {
    return this->densities;
}

glm::vec4* ParticleSystem::getTemperatures() {
    return this->temperatures;
}

glm::vec4* ParticleSystem::getForces() {
    return this->forces;
}

void ParticleSystem::setAccelerations(glm::vec4 *newAccelerations) {
    delete [] this->accelerations;
    this->accelerations = newAccelerations;
}

void ParticleSystem::setMasses(glm::vec4 *newMasses) {
    delete [] this->masses;
    this->masses = newMasses;
}

void ParticleSystem::setPositions(glm::vec4 *newPositions) {
    delete [] this->positions;
    this->positions = newPositions;
}

void ParticleSystem::setVelocities(glm::vec4 *newVelocities) {
    delete [] this->velocities;
    this->velocities = newVelocities;
}

void ParticleSystem::setTemperatures(glm::vec4 *newTemps) {
    delete [] this->temperatures;
    this->temperatures = newTemps;
}

std::ostream& operator<<(std::ostream& os, const ParticleSystem& system) {
    os << "Particle System with " << system.numParticles << " particles:" << std::endl;
    for (unsigned int i = 0; i < system.numParticles; ++i) {
        os << "Particle ID: " << i << std::endl;
        os << "Position: (" << system.positions[i].x << ", " << system.positions[i].y << ", " << system.positions[i].z << ")" << std::endl;
        os << "Velocity: (" << system.velocities[i].x << ", " << system.velocities[i].y << ", " << system.velocities[i].z << ")" << std::endl;
        os << "Acceleration: (" << system.accelerations[i].x << ", " << system.accelerations[i].y << ", " << system.accelerations[i].z << ")" << std::endl;
        os << "Mass: " << system.masses[i].x << std::endl;
    }
    return os;
}