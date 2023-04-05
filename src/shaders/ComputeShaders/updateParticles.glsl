#version 440 core

layout( local_size_x = 64, local_size_y =1, local_size_z = 1  ) in;

uniform int numParticles;
uniform float deltaTime;

layout(std430, binding=0) buffer positionsBuffer
{
    vec4 positions[];
};

layout(std430, binding=1) buffer velocitiesBuffer
{
    vec4 velocities[];
};

layout(std430, binding=2) buffer accelerationsBuffer
{
    vec4 accelerations[];
};

/**
 * Updates the particle position
 * Performs the leapfrog integration
 * @param index
 */
void updatePosition(uint index) {
    float dtDividedBy2 = deltaTime/2.f;

    // Compute velocity (i + 1/2)
    velocities[index] += accelerations[index] * dtDividedBy2;

    // Compute next position (i+1)
    positions[index] += velocities[index] * deltaTime;

    // Update acceleration (i+1)
    // TODO: new Acceleration
    accelerations[index] = vec4(0);

    // Compute next velocity (i + 1)
    velocities[index] += accelerations[index] * dtDividedBy2;
}



void main() {
    uint index = gl_GlobalInvocationID.x;

    if (index < numParticles) {
        updatePosition(index);
    }

}