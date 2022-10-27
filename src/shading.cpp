#include "texture.h"
#include "draw.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 position = ray.origin + ray.t * ray.direction;
    glm::vec3 lightDir = glm::normalize(lightPosition - position);
    float dot = glm::dot(glm::normalize(hitInfo.normal), lightDir);
    float eps = 1e-6;
    if (dot < eps) {
        return {0, 0, 0};
    }
    glm::vec3 view = ray.origin - position;
    glm::vec3 h = glm::normalize(lightDir + view);
    if (glm::dot(h, glm::normalize(hitInfo.normal)) < eps) {
        return {0, 0, 0};
    }
    // EXPERIMENT
    if (hitInfo.material.kdTexture && features.enableTextureMapping) {
        glm::vec3 texel = acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
        if (features.extra.enableBilinearTextureFiltering) {
            texel = bilinearInterpolation(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
        }
        return lightColor * hitInfo.material.ks * pow(glm::dot(glm::normalize(hitInfo.normal), h), hitInfo.material.shininess) + lightColor * texel * dot; 
    }
    return lightColor * hitInfo.material.ks * pow(glm::dot(glm::normalize(hitInfo.normal), h), hitInfo.material.shininess) + lightColor * hitInfo.material.kd * dot;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    glm::vec3 point = ray.origin + ray.direction * ray.t;
    //glm::vec3 r = glm::reflect(glm::normalize(ray.direction), glm::normalize(hitInfo.normal));
    glm::vec3 r = glm::normalize(ray.direction) - 2 * glm::dot(glm::normalize(hitInfo.normal), glm::normalize(ray.direction)) * glm::normalize(hitInfo.normal);
    Ray reflectionRay {point + float(1e-5)*r, r};
    //drawRay(Ray{point, hitInfo.normal, 1}, glm::vec3{0,1,0.5});
    return reflectionRay;
}