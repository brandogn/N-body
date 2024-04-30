// vertex shader
#version 440 core

layout(std430, binding=0) buffer positionsBuffer
{
    vec4 positions[];
};

layout(std430, binding=1) buffer velocitiesBuffer
{
    vec4 velocities[];
};

layout(std430, binding=3) buffer massesBuffer
{
    vec4 masses[];
};


out vec4 particleVelocity;

// 3d view stuff
uniform mat4 modelViewProjection;

uniform vec3 cameraPos;

uniform float worldSize;

uniform bool pointSize;

out float particleSize;
out float temp;

float getParticleSize(){
    // set the particle size based on the distance from the camera to the particle
    vec3 worldPos = positions[gl_VertexID].xyz;
    float particleSize = 4.0f / length(cameraPos - worldPos);

    float origin_dist_multiplier = abs(length(worldPos - vec3(2.5f))) * 3.0f;
    // float origin_dist_multiplier = clamp(
    //     abs(length(worldPos - vec3(2.5f))),
    //     0.0f,
    //     50.0f);

    // set the point size based on the particle size
    if(pointSize){
        return particleSize * worldSize * origin_dist_multiplier;
    }
    return origin_dist_multiplier;
}

float calcTemp() {
    float squaredVelocity = dot(velocities[gl_VertexID], velocities[gl_VertexID]);
    return ((2 * squaredVelocity) / (3 * 1.38e-23 * masses[gl_VertexID].x));
}


void main()
{
    // Set the point size
    gl_PointSize = getParticleSize();
    particleSize = gl_PointSize;

    // Set the position of the particle
    gl_Position = modelViewProjection * vec4(positions[gl_VertexID].xyz, 1.f);

    // Pass the velocity to the fragment shader
    particleVelocity = velocities[gl_VertexID];
    temp = calcTemp();
}