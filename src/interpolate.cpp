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
