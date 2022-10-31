#include "transparency.h"
#include "render.h"
#include "intersect.h"


glm::vec3 calculateColorTransparency(const Scene &scene, Ray &ray, const BvhInterface& bvh,  const Features &features, int rayDepth) {
    HitInfo hitInfo;
    Ray copyy = ray;
    bvh.intersect(ray, hitInfo, features);

    if (ray.t == std::numeric_limits<float>::max() || hitInfo.material.transparency == 1.0f) { 
        ray.t = std::numeric_limits<float>::max();
        return getFinalColor(scene, bvh, ray, features, rayDepth);
    }
    
    float opp = 1.0f - hitInfo.material.transparency;
    Ray newRay = {ray.origin + glm::vec3(ray.t + 1e-6) * ray.direction, ray.direction};
    ray.t = std::numeric_limits<float>::max();
    return hitInfo.material.transparency * getFinalColor(scene, bvh, ray, features, rayDepth) + opp * calculateColorTransparency(scene, newRay, bvh, features, rayDepth); 
}