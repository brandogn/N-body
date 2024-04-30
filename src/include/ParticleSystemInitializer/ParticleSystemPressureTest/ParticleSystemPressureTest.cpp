#include "ParticleSystemPressureTest.h"
#include <random>
#include <glm/gtc/random.hpp>
#include <chrono>

ParticleSystemPressureTest::ParticleSystemPressureTest(size_t numParticles) : totalParticles(numParticles){}

ParticleSystem* ParticleSystemPressureTest::generateParticles(glm::vec3 worldDimensions) {
    Particle* particles = new Particle[this->totalParticles];

    float length = glm::length(worldDimensions) / 20.f;

    // Determine radius of sphere based on world dimensions
    float radius = glm::length(worldDimensions) / 3.0f;

    // Determine particle mass based on world volume and total particles
    float worldVolume = worldDimensions.x * worldDimensions.y * worldDimensions.z;
    float particleMass = (worldVolume / this->totalParticles) * 1.0;

    
    float coreRatio = (1.f / 5.f);
    int numCoreParticles = totalParticles * coreRatio;

    float shockRadius = radius * coreRatio; // Radius from center to differentiate between the inner and outer layers
    float shockStrength = 200.0f; // Initial outward velocity for outer particles

    for (int i = 0; i < numCoreParticles; ++i) {
        glm::vec3 position = glm::ballRand(radius * coreRatio);
        float distToCenter = glm::length(position);
        glm::vec3 initialVel = (distToCenter > shockRadius) ? shockStrength * glm::normalize(position) : glm::vec3(0.0f);

        // Shift particle position by half of the world dimensions
        position += 0.5f * worldDimensions;
        //(worldDimensions.z == 0) ? glm::linearRand(glm::vec3(-length, -length, 0.f), glm::vec3(length, length, 0.f)) : glm::linearRand(glm::vec3(-length, -length, -length), glm::vec3(length, length, length));
        particles[i] = Particle(position, initialVel, particleMass * 1000.0f);
    }

    for (int i = numCoreParticles; i < totalParticles; ++i) {
        glm::vec3 position = glm::ballRand(radius);
        float distToCenter = glm::length(position);
        glm::vec3 initialVel = (distToCenter > shockRadius) ? shockStrength * glm::normalize(position) : glm::vec3(0.0f);

        // Shift particle position by half of the world dimensions
        position += 0.5f * worldDimensions;
        //(worldDimensions.z == 0) ? glm::linearRand(glm::vec3(-length, -length, 0.f), glm::vec3(length, length, 0.f)) : glm::linearRand(glm::vec3(-length, -length, -length), glm::vec3(length, length, length));
        particles[i] = Particle(position, initialVel, particleMass);
    }

    // for (int i = 0; i < totalParticles; ++i) {
    //     particles[i] = Particle(glm::vec3(0.0f), glm::vec3(0.0f), 0.0f);
    // }

    std::vector<Particle> particleVector(particles, particles + this->totalParticles);

    delete[] particles; // free the allocated memory

    return new ParticleSystem(particleVector);
}