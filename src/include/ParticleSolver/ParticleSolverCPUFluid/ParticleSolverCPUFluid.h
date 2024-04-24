
#include "ParticleSolver.h"

#ifndef N_BODY_PARTICLESOLVERCPUFLUID_H
#define N_BODY_PARTICLESOLVERCPUFLUID_H

class ParticleSolverCPUFluid : public ParticleSolver {
public:
  ParticleSolverCPUFluid(float timeStep, float squaredSoftening);
  void updateParticlePositions(ParticleSystem *particles) override;
  bool usesGPU() override;
  float getSquaredSoftening() override;

protected:
  float squaredSoftening;
  float G;
  float timeStep;
  // density_map
  void computeDensityMap(ParticleSystem *particles);
  void computeGravityForce(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeFluidForce(ParticleSystem *particles,
                          const unsigned int particleId);
};

#endif // N_BODY_PARTICLESOLVERCPUFLUID_H
