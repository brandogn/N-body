#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "ParticleSystemFile.h"

ParticleSystemFile::ParticleSystemFile(std::string filePath) : path(filePath){}

ParticleSystem* ParticleSystemFile::generateParticles(glm::vec3 worldDimensions) {
    std::ifstream file(path);
    if (!file) {
        std::cerr << "Failed to open file: " << path << std::endl;
        exit(EXIT_FAILURE); // Terminate the program
    }

    std::string line;
    std::getline(file, line);
    int totalParticles = 0;
    if (line.find("Particle System with ") != std::string::npos) {
        // Extract total particles from the first line
        std::string particlesCountStr = line.substr(line.find("with ") + 5);
        totalParticles = std::stoi(particlesCountStr);
    } else {
        std::cerr << "Invalid file format. Unable to determine total particles." << std::endl;
        exit(EXIT_FAILURE); // Terminate the program
    }

    Particle* particles = new Particle[totalParticles];

    int currentParticleIndex = 0;
    glm::vec3 position, velocity, acceleration;
    float mass;

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        if (line.find("Particle ID") != std::string::npos) {
            // Skip particle ID line
            continue;
        } else if (line.find("Position:") != std::string::npos) {
            // Extract particle position
            iss.ignore(std::numeric_limits<std::streamsize>::max(), '(');
            std::string posX, posY, posZ;
            std::getline(iss, posX, ',');
            std::getline(iss, posY, ',');
            std::getline(iss, posZ, ')');
            position.x = std::stof(posX);
            position.y = std::stof(posY);
            position.z = std::stof(posZ);
            std::cout << "Position: " << position.x << ", " << position.y << ", " << position.z << std::endl;
        } else if (line.find("Velocity:") != std::string::npos) {
            // Extract particle velocity
            iss.ignore(std::numeric_limits<std::streamsize>::max(), '(');
            std::string velX, velY, velZ;
            std::getline(iss, velX, ',');
            std::getline(iss, velY, ',');
            std::getline(iss, velZ, ')');
            velocity.x = std::stof(velX);
            velocity.y = std::stof(velY);
            velocity.z = std::stof(velZ);
            std::cout << "Velocity: " << velocity.x << ", " << velocity.y << ", " << velocity.z << std::endl;
        }
        else if (line.find("Acceleration:") != std::string::npos) {
            // Extract particle acceleration
            iss.ignore(std::numeric_limits<std::streamsize>::max(), '(');
            std::string accX, accY, accZ;
            std::getline(iss, accX, ',');
            std::getline(iss, accY, ',');
            std::getline(iss, accZ, ')');
            acceleration.x = std::stof(accX);
            acceleration.y = std::stof(accY);
            acceleration.z = std::stof(accZ);
            std::cout << "Acceleration: " << acceleration.x << ", " << acceleration.y << ", " << acceleration.z << std::endl;
        } else if (line.find("Mass:") != std::string::npos) {
            // Extract particle mass
            std::string massStr = line.substr(line.find(":") + 1);
            mass = std::stof(massStr);
            std::cout << "Mass: " << mass << std::endl;

            // Create Particle instance with extracted data
            particles[currentParticleIndex] = Particle(position, velocity, acceleration, mass);
            currentParticleIndex++;

            if (currentParticleIndex >= totalParticles) {
                break; // Stop reading the file if all particles have been processed
            }
        }
    }
    file.close();

    std::vector<Particle> particleVector(particles, particles + totalParticles);

    delete[] particles; // free the allocated memory

    return new ParticleSystem(particleVector);
}