#pragma once


#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include "common.h"
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <random>
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include "common.h"
#include "bvh_interface.h"

/*

    EXTRA FEATURE: DEPTH OF FIELD

*/

glm::vec3 getFocalPoint(Ray ray, float focalLength);
glm::vec3 generateNewOrigin(glm::vec3 origin, float aperture);
std::vector<glm::vec3> generateSamples(glm::vec3 origin, float aperture, int samples, std::default_random_engine &rng, std::uniform_real_distribution<float> &dist);
std::vector<Ray> sampledRays(glm::vec3 point, std::vector<glm::vec3> origins); 
void debugDOF(Ray ray, std::vector <Ray>& defaultRaysDOF, glm::vec3& Lo, int a, int &b);
glm::vec3 debugDepthOfField(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth);