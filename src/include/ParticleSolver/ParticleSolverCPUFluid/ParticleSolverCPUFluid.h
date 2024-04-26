
#include "ParticleSolver.h"
#include "GridCPU.h"

#ifndef N_BODY_PARTICLESOLVERCPUFLUID_H
#define N_BODY_PARTICLESOLVERCPUFLUID_H

#define PI 3.14159265358979323846


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

  float densityKernel(float distance, float radius);
  float densityDerivative(float distance, float radius);
  float nearDensityKernel(float distance, float radius);
  float nearDensityDerivative(float distance, float radius);
  float pressureFromDensity(float density);
  float nearPressureFromDensity(float density);

  void computeDensityMap(ParticleSystem *particles,
                           const unsigned int particleId);
  void computePressureForce(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeGravityForce(ParticleSystem *particles,
                           const unsigned int particleId);
  void computeFluidForce(ParticleSystem *particles,
                          const unsigned int particleId);
};

#endif // N_BODY_PARTICLESOLVERCPUFLUID_H
