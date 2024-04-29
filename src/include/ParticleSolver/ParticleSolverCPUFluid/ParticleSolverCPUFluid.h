
#include "ParticleSolver.h"
#include "GridCPU.h"
#include "FluidMath.h"

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
  float targetDensity;
  float pressureMultiplier;
  float nearPressureMultiplier;
  float viscosityStrength;

  float pressureFromDensity(float density);
  float nearPressureFromDensity(float density);

  void computeDensityMap(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeTemperatures(ParticleSystem *particles,
                           const unsigned int particleId);
  void computePressureForce(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeViscosityForce(ParticleSystem *particles, 
                           const unsigned int particleId);
  void computeGravityForce(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeFluidForce(ParticleSystem *particles,
                          const unsigned int particleId);
};

#endif // N_BODY_PARTICLESOLVERCPUFLUID_H
