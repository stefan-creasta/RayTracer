#include "texture.h"
#include "draw.h"
#include <cmath>
#include <glm/geometric.hpp>
#include "shading.h"
#include <iostream>

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
        //ImageMipMap mipmap = getMipMap(*hitInfo.material.kdTexture.get());
        if (features.extra.enableBilinearTextureFiltering) {

            texel = bilinearInterpolation(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
            glm::vec2 debugUV = getUVForBilinear(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
            glm::vec3 w = glm::normalize(hitInfo.normal);
            glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
            glm::vec3 u = glm::normalize(glm::cross(t, w));
            glm::vec3 v = glm::normalize(glm::cross(w, u));
            Ray rayU = Ray { position, u, debugUV.x / 5.0f };
            Ray rayV = Ray { position, v, debugUV.y / 5.0f };
            Ray rayOU = Ray { position, -u, (1 - debugUV.x) / 5.0f };
            Ray rayOV = Ray { position, -v, (1 - debugUV.y) / 5.0f };
            drawRay(rayU, glm::vec3 { 1.0f, 0.5f, 0.0f });
            drawRay(rayV, glm::vec3 { 1.0f, 0.0f, 0.5f });
            drawRay(rayOU, glm::vec3 { 1.0f, 0.5f, 0.0f });
            drawRay(rayOV, glm::vec3 { 1.0f, 0.0f, 0.5f });
            //std::cout << rayU.origin.x << std::endl;
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

/** glm::vec3 trilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features, const Ray& ray, const BvhInterface& bvh)
{
    // Get the mipmap first
    ImageMipMap mipmap = getMipMap(image);
    glm::vec3 point = ray.origin + ray.t * ray.direction;
    Ray rayCopy = ray;
    HitInfo hitInfo;
    // Get the triangle which was hit, in order to calculate the texCoord and 3D coordinates for a corner
    triangleMeshPair getT = bvh.getTriangleIntersection(rayCopy, hitInfo, features);
    Mesh mesh = getT.mesh;
    glm::uvec3 tri = getT.tri;
    // If a sphere was hit first, we compute bilinear interpolation
    if (tri == glm::uvec3(-1000000)) {
        return bilinearInterpolation(image, texCoord, features);
    }
    auto v = mesh.vertices[tri[0]];
    if (point == v.position) {
        v = mesh.vertices[tri[1]];
    }
    // Computing the derivative
    // First we calculate the 2D distance between the texture coordinates
    glm::vec2 texCoord2 = v.texCoord;
    float distanceTexX = texCoord.x - texCoord2.x;
    float distanceTexY = texCoord.y - texCoord2.y;
    // Then we approximate the screen space
    Ray newRay = ray;
    newRay.direction = glm::normalize(v.position - ray.origin) * glm::length(ray.direction);
    newRay.t = ray.t;
    glm::vec3 newPoint = newRay.t * newRay.direction + newRay.origin;
    float distanceScreenX = point.x - newPoint.x;
    float distanceScreenY = point.y - newPoint.y;
    // The derivative is the distanceScreen / distanceTex
    float derivativeX = 1.0, derivativeY = 1.0;
    // If the distance for coordinates is zero, the derivative will converge to 1, so it stays 1
    if (distanceTexX != 0.0f) {
        derivativeX = distanceScreenX / distanceTexX;
    }
    // Same for y axis
    if (distanceTexY != 0.0f) {
        derivativeY = distanceScreenY / distanceTexY;
    }
    // Taking the maximum derivative
    float maxDerivative = derivativeX;
    if (maxDerivative < derivativeY) {
        maxDerivative = derivativeY;
    }
    // Finding the log
    float k = std::log2(maxDerivative);
    float k0 = std::floor(k);
    float k1 = k0 + 1;
    float a = k1 - k;
    // The level cannot be negative, so we approximate the trilinear interpolation to a bilinear interpolation for the first level in the mipmap (level 0)
    if (k0 < 0) {
        return bilinearInterpolation(image, texCoord, features);
    }
    glm::vec3 c0 = bilinearInterpolationForMipMap(mipmap, k0, texCoord, features);
    glm::vec3 c1 = bilinearInterpolationForMipMap(mipmap, k1, texCoord, features);
    return a * c0 + (1 - a) * c1;
}**/