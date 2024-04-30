#include <iostream>
#include <cmath>

#define PI 3.14159265358979323846
#define BOLTZMANN 1.38e-23

class FluidMath {
    public:
        static float densityKernel(float distance, float radius) {
            if (distance < radius) {
                float scale = 15 / (2 * PI * std::pow(std::abs(radius), 5));
                float v = radius - distance;
                return v * v * scale;
            }
            return 0;
        }
        static float densityDerivative(float distance, float radius) {
            if (distance <= radius) {
                float scale = 15 / (PI * std::pow(std::abs(radius), 5));
                float v = radius - distance;
                return -v * scale;
            }
            return 0;
        }
        static float nearDensityKernel(float distance, float radius) {
            if (distance < radius) {
                float scale = 15 / (PI * std::pow(std::abs(radius), 6));
                float v = radius - distance;
                return v * v * v * scale;
            }
            return 0;
        }
        static float nearDensityDerivative(float distance, float radius) {
            if (distance <= radius) {
                float scale = 45 / (PI * std::pow(std::abs(radius), 6));
                float v = radius - distance;
                return -v * v * scale;
            }
            return 0;
        }
        static float viscosityKernel(float distance, float radius) {
            if (distance < radius) {
                float scale = 315 / (64 * PI * std::pow(std::abs(radius), 9));
                float v = radius * radius - distance * distance;
                return v * v * v;
            }
            return 0;
        }
        static float computeTemperature(float mass, float velocity) {
            return ((2 * velocity * velocity) / (3 * BOLTZMANN * mass));
        }
};