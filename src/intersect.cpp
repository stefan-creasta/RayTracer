#include "intersect.h"
#include "interpolate.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>
#include <iostream>

constexpr float maxFloat = std::numeric_limits<float>::infinity();

bool isDegenerateTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    return glm::length(glm::cross(v2 - v0, v1 - v0)) == 0.f;
}

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    if (isDegenerateTriangle(v0, v1, v2))
        return false;

    const float alpha = glm::dot(n, glm::cross(v2 - v1, p - v1));
    const float beta = glm::dot(n, glm::cross(p - v0, v2 - v0));
    const float gamma = glm::dot(n, glm::cross(v1 - v0, p - v0));

    return !(alpha < 0 || beta < 0 || gamma < 0);
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    if (glm::dot(plane.normal, ray.direction) == 0.f)
        return false;
    const float rayDistanceToOrigin = glm::dot(plane.normal, ray.origin);
    const float t = (plane.D - rayDistanceToOrigin) / glm::dot(plane.normal, ray.direction);
    if (t <= 0.f) {
        return false;
    } else {
        if (t < ray.t) {
            ray.t = t;
            return true;
        } else {
            return false;
        }
    }
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

    float D = glm::dot(v0, normal);

    Plane plane { D, normal };
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    const float prevt = ray.t;

    if (isDegenerateTriangle(v0, v1, v2))
        return false;

    const Plane plane = trianglePlane(v0, v1, v2);

    if (intersectRayWithPlane(plane, ray)) {
        const glm::vec3 n = plane.normal;
        const glm::vec3 p = ray.origin + ray.t * ray.direction;

        const float alpha = glm::dot(n, glm::cross(v2 - v1, p - v1));
        const float beta = glm::dot(n, glm::cross(p - v0, v2 - v0));
        const float gamma = glm::dot(n, glm::cross(v1 - v0, p - v0));

        if (!(alpha < 0 || beta < 0 || gamma < 0)) {
            const float invTriangleArea = 1 / glm::dot(n, glm::cross(v1 - v0, v2 - v0));
            hitInfo.barycentricCoord = glm::vec3 { alpha * invTriangleArea, beta * invTriangleArea, gamma * invTriangleArea };
            hitInfo.normal = n;
            return true;
        }
    }

    ray.t = prevt;
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float a = glm::dot(ray.direction, ray.direction);
    float b = 2 * glm::dot(ray.direction, (ray.origin - sphere.center));
    float c = glm::dot(sphere.center, sphere.center) + glm::dot(ray.origin, ray.origin) - 2 * glm::dot(ray.origin, sphere.center) - sphere.radius * sphere.radius;
    float delta = b * b - 4 * a * c;
    if (delta < 0) {
        return false;
    }
    float t1 = (-b + sqrt(delta)) / (2 * a);
    float t2 = (-b - sqrt(delta)) / (2 * a);
    glm::vec3 v1 = ray.origin + t1 * ray.direction;
    glm::vec3 v2 = ray.origin + t2 * ray.direction;

    if (t2 < 0) {
        if (ray.t > t2) {
            ray.t = std::min(t1, ray.t);
        } else {
            return false;
        }
    } else {
        if (t1 > 0.f && ray.t > t1) {
            ray.t = std::min(t2, ray.t);
        } else {
            return false;
        }
    }
    glm::vec3 hit = ray.origin + ray.t * ray.direction;
    hitInfo.normal = glm::normalize(hit - sphere.center);
    hitInfo.material = sphere.material;
    return true;
}


/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float first = 0.0f, second = maxFloat;
    for (int i = 0; i < 3; i++) {
        float inverseDirection = 1.0f / ray.direction[i];
        float t0 = (box.lower[i] - ray.origin[i]) * inverseDirection;
        float t1 = (box.upper[i] - ray.origin[i]) * inverseDirection;

        if (t0 > t1)
            std::swap(t0, t1);

        if (t0 > second || t1 < first)
            return false;

        if (t0 > first)
            first = t0;
        if (t1 < second)
            second = t1;
    }

    if (first <= 0.0f)
        first = second;

    if (ray.t < first) { 
        return false;
    } else {
        ray.t = first;
        return true;
    }
}
