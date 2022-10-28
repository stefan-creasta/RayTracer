#include "dof.h"
#include "draw.h"
#include "render.h"
#include <iostream>

/*
    Focal point where all the rays should intersect.
*/
glm::vec3 getFocalPoint(Ray ray, float focalLength) {
    return ray.origin + focalLength * ray.direction; 
}


/*
    Generate new random origin position.
*/
glm::vec3 generateNewOrigin(glm::vec3 origin, float aperture) {
    float xshift = -aperture / 2 + static_cast<float>(rand()) / (RAND_MAX / aperture);
    float yshift = -aperture / 2 + static_cast<float>(rand()) / (RAND_MAX / aperture);
    float zshift = -aperture / 2 + static_cast<float>(rand()) / (RAND_MAX / aperture);
    return origin + glm::vec3 {xshift, yshift, zshift}; 
}

/*
    Generate multiple samples of origin with some aperture.
*/
std::vector<glm::vec3> generateSamples(glm::vec3 origin, float aperture, int samples, std::default_random_engine &rng, std::uniform_real_distribution<float> &dist) {

    std::vector<glm::vec3> origins;
    for (int i = 0; i < samples; i++) {
    
       float xshift = dist(rng);
       float yshift = dist(rng);
       float zshift = dist(rng);
       origins.push_back(origin + glm::vec3{xshift, yshift, zshift});
    }
    return origins;
}

/*
    Generate multiple rays for depth of field effect.
*/
std::vector<Ray> sampledRays(glm::vec3 point, std::vector<glm::vec3> origins) {
    std::vector<Ray> rays;
    for (auto origin : origins) {
        rays.push_back(Ray {origin, glm::normalize((point - origin))});
    }
    return rays;
}

glm::vec3 debugDepthOfField(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth) {

    float aperture = 0.2;
    float focalLength = 2;
    std::default_random_engine rng{42};
    std::uniform_real_distribution<float> dist(-aperture/2, aperture/2);
    glm::vec3 col = glm::vec3 {0, 0, 0};
    int samples = 100;
    
    glm::vec3 focalPoint = getFocalPoint(ray, focalLength);
    std::vector <Ray> rays = sampledRays(focalPoint, generateSamples(ray.origin, aperture, samples, rng, dist));
    std::cout << rayDepth << "\n";
    for (auto r : rays) {
        col += getFinalColor(scene, bvh, r, features, rayDepth);
    }

    col = col * glm::vec3(1.0f / (rays.size())); 
    return col;
}