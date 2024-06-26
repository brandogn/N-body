#include <glm/glm.hpp>
#include <vector>
#include <ostream>
#include <Particle.h>

#ifndef N_BODY_PARTICLESYSTEM_H
#define N_BODY_PARTICLESYSTEM_H


class ParticleSystem {
public:
    ParticleSystem(std::vector<Particle> &particles);
    ParticleSystem(ParticleSystem* other);
    void updateParticlePosition(unsigned int particleId, float deltaTime);
    unsigned int size() const;
    glm::vec4* getPositions();
    glm::vec4* getVelocities();
    glm::vec4* getAccelerations();
    glm::vec4* getMasses();
    glm::vec4* getDensities();
    glm::vec4* getTemperatures();
    glm::vec4* getForces();
    void setMasses(glm::vec4* newMasses);
    void setPositions(glm::vec4* newPositions);
    void setAccelerations(glm::vec4* newAccelerations);
    void setVelocities(glm::vec4* newVelocities);
    void setTemperatures(glm::vec4* newTemps);
    friend std::ostream& operator<<(std::ostream& os, const ParticleSystem& system);

protected:

    unsigned int numParticles;
    glm::vec4* positions;
    glm::vec4* accelerations;
    glm::vec4* velocities;
    glm::vec4* masses;
    glm::vec4* densities;
    glm::vec4* temperatures;
    glm::vec4* forces;

    // Declare ParticleSimulation as a friend class
    friend class ParticleSimulation;

};


#endif //N_BODY_PARTICLESYSTEM_H
