//fragshader
#version 440 core
layout (location = 0) out vec4 FragColor;
layout (location = 1) out vec4 BrightColor;

uniform float worldSize;

in vec4 particleVelocity;

in float particleSize;
in float temp;

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec4 getParticleColor(){
    // Define a maximum velocity value
    float maxVelocity = worldSize*1.5;

    // Compute the magnitude of the particle's velocity
    float velocityMagnitude = length(particleVelocity);

    // Compute a normalized velocity value between 0 and 1
    float normalizedVelocity = clamp(velocityMagnitude / maxVelocity, 0.0, 1.0);

    // // Define three colors for the gradient (e.g. red, orange, yellow)
    // vec3 colorLow = vec3(0.0, 0.0, 1.0); // blue (slowest)
    // vec3 colorMid = vec3(1.0, 0.5, 0.0); // orange (middle)
    // vec3 colorHigh = vec3(1.0, 1.0, 0.0); // yellow (fastest)

    // // Interpolate between the three colors based on the normalized velocity value
    // float smoothNormalizedVelocity1 = smoothstep(0.0, 0.5, normalizedVelocity);
    // float smoothNormalizedVelocity2 = smoothstep(0.5, 1.0, normalizedVelocity);

    // vec3 color = mix(colorLow, colorMid, smoothNormalizedVelocity1);
    // color = mix(color, colorHigh, smoothNormalizedVelocity2);

    // return vec4(color, 1); // pass the velocity to the fragment shader


    
    // Define a maximum temp value
    float maxtemp = 100000000.0f;

    // Compute a normalized temp value between 0 and 1
    float normalizedtemp = clamp(temp / maxtemp, 0.0, 1.0);

    vec3 hsv = vec3(normalizedtemp, 1.0f, normalizedVelocity);


    return vec4(hsv2rgb(hsv), 1.0f);
}


vec4 getParticleBrightColor(){
    // Glow the particle
    return vec4(FragColor.rgb, 1.0);
}

float getAlpha(float dist) {
    return (0.25f-abs(dist)) * 4.0f / particleSize;
}

void main()
{
    // Make GL_POINTS circular
    vec2 pos = gl_PointCoord.xy-0.5;
    float dist = dot(pos,pos);
    if(dist > 0.25 && particleSize > 3.0){
        discard;
    }

    // Set the color of the particle
    FragColor = getParticleColor();

    // Render in the second framebuffer the bright particles
    BrightColor = getParticleBrightColor();

    float a;
    if (dist < 0.00005f) {
        a = 1.0f;
    } else {
        a = getAlpha(dist);
    }

    FragColor[3] = a;
    BrightColor[3] = a;
    // FragColor   = vec4(a);
    // BrightColor = vec4(a);
}