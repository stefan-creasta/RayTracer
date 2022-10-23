#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.
    float tArea = glm::length(glm::cross(v1 - v2, v2 - v0));
    float a = glm::length(glm::cross(v1 - p, v2 - p)) / tArea;
    float b = glm::length(glm::cross(v0 - p, v2 - p)) / tArea;
    float g = glm::length(glm::cross(v1 - p, v0 - p)) / tArea;
    return glm::vec3{a, b, g};
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    // TODO: implement this function.
    return barycentricCoord[0] * n0 + barycentricCoord[1] * n1 + barycentricCoord[2] * n2;
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord[0] * t0 + barycentricCoord[1] * t1 + barycentricCoord[2] * t2;
}

void interpolateNormalDebug(const Vertex v0, const Vertex v1, const Vertex v2, const Ray ray, const HitInfo hitInfo) {
    drawRay(Ray{v0.position, v0.normal, 1}, glm::vec3{1, 0.2, 0.4});
    drawRay(Ray{v1.position, v1.normal, 1}, glm::vec3{0.5, 1, 0});
    drawRay(Ray{v2.position, v2.normal, 1}, glm::vec3{0, 0.9, 1});
    drawRay(Ray{ray.origin + ray.t * ray.direction, hitInfo.normal, 1}, glm::vec3{0, 1, 0}); 
}