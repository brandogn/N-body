
#include "ParticleSystemInitializer.h"

#ifndef PARTICLESYSTEMPRESSURETEST_H
#define PARTICLESYSTEMPRESSURETEST_H


class ParticleSystemPressureTest: public ParticleSystemInitializer{

public:
    ParticleSystemPressureTest(size_t numParticles);
    ParticleSystem* generateParticles(glm::vec3 worldDimensions);

private:
    size_t totalParticles;
};


#endif //PARTICLESYSTEMPRESSURETEST_H
