
#include "ParticleSolver.h"
#include "GridCPU.h"

#ifndef N_BODY_PARTICLESOLVERCPUFLUID_H
#define N_BODY_PARTICLESOLVERCPUFLUID_H

class ParticleSolverCPUFluid : public ParticleSolver {
public:
  ParticleSolverCPUFluid(GridCPU *grid, float timeStep, float squaredSoftening);
  ~ParticleSolverCPUFluid();
  void updateParticlePositions(ParticleSystem *particles) override;
  bool usesGPU() override;
  float getSquaredSoftening() override;

protected:
  float squaredSoftening;
  float G;
  float timeStep;
  GridCPU *grid;
  // density_map
  float smoothingRadius;
  float smoothingKernelDensity(float distance, float radius);
  void computeDensityMap(ParticleSystem *particles);
  void computeGravityForce(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeFluidForce(ParticleSystem *particles,
                          const unsigned int particleId);
};

#endif // N_BODY_PARTICLESOLVERCPUFLUID_H
