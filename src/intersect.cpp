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


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    glm::vec3 bary = computeBarycentricCoord(v0, v1, v2, p);
    float a = bary[0];
    float b = bary[1];
    float g = bary[2];
    
    if (a < 0.0f || b < 0.0f || g < 0.0f) {
        return false;
    }

    float eps = 1e-6;
    return (a + b + g <= 1.0f + eps && a + b + g >= 1.0f - eps);
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float check = glm::dot(plane.normal, glm::normalize(ray.origin + ray.direction));
    float eps = 1e-6;
    if (check <= eps && check >= -eps) {
        return false;
    }
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / (glm::dot(ray.direction, plane.normal));
    if (t < 0) {
        return false;
    }
    if (t > ray.t) {
        return false;
    }
    ray.t = t;
    return true;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 l1 = v1 - v0, l2 = v2 - v0;
    plane.normal = glm::cross(l1, l2);
    if (plane.normal == glm::vec3 {0, 0, 0}) {
        plane.normal = glm::vec3 {1, 1, 1};
    }
    plane.normal = glm::normalize(plane.normal);
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{

   //return false; // CAREFULLL
   Plane plane = trianglePlane(v0, v1, v2);
   float temp = ray.t;
   if (!intersectRayWithPlane(plane, ray)) {
        return false;
   }
   if (!pointInTriangle(v0, v1, v2, plane.normal, ray.origin + ray.t * ray.direction)) {
        ray.t = temp;
        return false;
   }
   ray.t = std::min(ray.t, temp);
   glm::vec3 hit = ray.origin + ray.t * ray.direction;
   hitInfo.normal = plane.normal;
   hitInfo.barycentricCoord = computeBarycentricCoord(v0, v1, v2, hit); 
   return true;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float a = glm::dot(ray.direction, ray.direction);
    float b = 2 * glm::dot(ray.direction, (ray.origin - sphere.center));
    float c = glm::dot(sphere.center, sphere.center) + glm::dot(ray.origin, ray.origin) - 2 * glm::dot(ray.origin, sphere.center) - sphere.radius * sphere.radius;
    float delta = b*b - 4*a*c;
    if (delta < 0) {
        return false;
    }
    float t1 = (-b + sqrt(delta)) / (2*a);
    float t2 = (-b - sqrt(delta)) / (2*a);
    glm::vec3 v1 = ray.origin + t1 * ray.direction;
    glm::vec3 v2 = ray.origin + t2 * ray.direction;

    if (t2 < 0) {
        ray.t = std::min(t1, ray.t);
    }
    else {
        ray.t = std::min(t2, ray.t);
    }
    glm::vec3 hit = ray.origin + ray.t * ray.direction;
    hitInfo.normal = glm::normalize(hit - sphere.center); 
    return true; 
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float first = 0.0f, second = MAXFLOAT;
    for (int i = 0; i < 3; i++) {
        float inverseDirection = 1.0f / ray.direction[i];
        float t0 = (box.lower[i] - ray.origin[i]) * inverseDirection;
        float t1 = (box.upper[i] - ray.origin[i]) * inverseDirection;

        if (t0 > t1) std::swap(t0, t1);

        if (t0 > second || t1 < first) return false;

        if (t0 > first) first = t0;
        if (t1 < second) second = t1;
    }

    if (first <= 0.0f) first = second; 
    ray.t = std::min(ray.t, first);
    return true;
}
