#include "../../lib/glad/glad.h"
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include "../include/ParticleSystem/AbstractClass/ParticleSystem.cpp"
#include "../include/OpenGLRenderer/OpenGLRenderer.cpp"
#include "../include/ParticleSystem/ParticleSystemCPU/ParticleSystemCPU.cpp"
#include "../include/ParticleSystem/ParticleSystemGPU/ParticleSystemGPU.cpp"
#include "../include/enums/enums.h"
#include "../include/ArgumentsParser/ArgumentsParser.cpp"


int main(int argc, char *argv[])
{
    // Get the arguments
    ArgumentsParser args(argc, argv);

    OpenGLRenderer renderer(600, 600, "N-body simulation", false, true);

    ParticleSystem* particleSystem;
    switch (args.getVersion()){
        case Version::PP_CPU_SEQUENTIAL:
            particleSystem = new ParticleSystemCPU(args.getNumParticles());
            break;
        case Version::PP_CPU_PARALLEL:
            std::cerr << "Not yet implemented \n";
            exit(EXIT_FAILURE);

        case Version::PP_GPU_PARALLEL:
            particleSystem = new ParticleSystemGPU(args.getNumParticles());
            break;
        default:
            exit(EXIT_FAILURE);
    }

    renderer.render_loop(particleSystem);

    delete particleSystem;




}
