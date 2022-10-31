#pragma once
#include "common.h"
#include <framework/ray.h>
extern float degreeBlur;
extern int numberOfRays;

struct ImageMipMap {
    std::vector<int> width;
    std::vector<int> height;
    std::vector<std::vector<glm::vec3>> pixels;
};

// Compute the shading at the intersection point using the Phong model.
const glm::vec3 computeShading (const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo);

// Given a ray and a normal (in hitInfo), compute the reflected ray in the specular direction (mirror direction).
const Ray computeReflectionRay (Ray ray, HitInfo hitInfo);

ImageMipMap getMipMap(const Image& image);
// Given a ray and a normal (in hitInfo), compute many reflected rays to compute glossy effect.
std::vector<Ray> glossyRays(Ray reflection);
