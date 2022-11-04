#include "texture.h"
#include "draw.h"
#include <cmath>
#include <glm/geometric.hpp>
#include "shading.h"
#include <random>
#include <iostream>
#include <glm/gtx/common.hpp>
int numberOfRays = 5;
bool showMipmapLevel = false;
int mipmapLevel = 0;

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 position = ray.origin + ray.t * ray.direction;
    glm::vec3 lightDir = glm::normalize(lightPosition - position);
    float dot = glm::dot(glm::normalize(hitInfo.normal), lightDir);
    float eps = 1e-6;
    if (dot < 0.0f) {
        return {0, 0, 0};
    }
    // EXPERIMENT
    auto reflectRay = computeReflectionRay({lightPosition, lightDir}, hitInfo);
    auto d = glm::abs(glm::dot(glm::normalize(reflectRay.direction), glm::normalize(ray.direction)));
    auto specular = lightColor * hitInfo.material.ks * glm::pow(d, hitInfo.material.shininess);
    if (hitInfo.material.kdTexture && features.enableTextureMapping) {
        glm::vec3 texel = acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
        //ImageMipMap mipmap = getMipMap(*hitInfo.material.kdTexture.get());
        if (features.extra.enableBilinearTextureFiltering) {
            if (features.extra.enableMipmapTextureFiltering) {
                texel = trilinearInterpolation(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features, ray, hitInfo);
            } else {

                texel = bilinearInterpolation(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
                glm::vec2 debugUV = getUVForBilinear(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
                glm::vec3 w = glm::normalize(hitInfo.normal);
                glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
                glm::vec3 u = glm::normalize(glm::cross(t, w));
                glm::vec3 v = glm::normalize(glm::cross(w, u));
                Ray rayU = Ray { position, -v, debugUV.x / 5.0f };
                Ray rayV = Ray { position, -u, debugUV.y / 5.0f };
                Ray rayOU = Ray { position, v, (1 - debugUV.x) / 5.0f };
                Ray rayOV = Ray { position, u, (1 - debugUV.y) / 5.0f };
                drawRay(rayU, glm::vec3 { 1.0f, 0.5f, 0.0f });
                drawRay(rayV, glm::vec3 { 1.0f, 0.0f, 0.5f });
                drawRay(rayOU, glm::vec3 { 1.0f, 0.5f, 0.0f });
                drawRay(rayOV, glm::vec3 { 1.0f, 0.0f, 0.5f });
                // std::cout << rayU.origin.x << std::endl;
            }
        }
        return lightColor * texel * dot + specular;
    }
    return lightColor * hitInfo.material.kd * dot + specular;
}

float getRandomVal2()
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float randomVal = distribution(generator);
    float randomVal2 = distribution(generator);
    // std::cout << randomVal << randomVal2 << std::endl;
    return randomVal;
}

Ray returnGlossyRay(Ray reflection, float degreeBlur)
{
    glm::vec3 w = glm::normalize(reflection.direction);
    glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
    glm::vec3 u = glm::normalize(glm::cross(t, w));
    glm::vec3 v = glm::normalize(glm::cross(w, u));
    float ua = -degreeBlur / 2.0f + degreeBlur * getRandomVal2();
    float va = -degreeBlur / 2.0f + degreeBlur * getRandomVal2();
    Ray rr = reflection;
    rr.direction += ua * u + va * v;
    glm::vec3 p1 = reflection.origin + u * degreeBlur / 10.0f + v * degreeBlur / 10.0f;
    glm::vec3 p2 = reflection.origin + u * degreeBlur / 10.0f - v * degreeBlur / 10.0f;
    glm::vec3 p3 = reflection.origin - u * degreeBlur / 10.0f + v * degreeBlur / 10.0f;
    glm::vec3 p4 = reflection.origin - u * degreeBlur / 10.0f - v * degreeBlur / 10.0f;
    Vertex v1 = Vertex(p1, reflection.direction);
    Vertex v2 = Vertex(p2, reflection.direction);
    Vertex v3 = Vertex(p3, reflection.direction);
    Vertex v4 = Vertex(p4, reflection.direction);
    drawTriangle(v1, v2, v3);
    drawTriangle(v2, v3, v4);
    return rr;
}

const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    glm::vec3 point = ray.origin + ray.direction * ray.t;
    glm::vec3 r = glm::normalize(ray.direction) - 2 * glm::dot(glm::normalize(hitInfo.normal), glm::normalize(ray.direction)) * glm::normalize(hitInfo.normal);
    Ray reflectionRay {point + float(1e-5)*r, r};
    return reflectionRay;
}

std::vector<Ray> glossyRays(Ray reflection, float degreeBlur)
{
    std::vector<Ray> rays;
    for (int i = 1; i <= numberOfRays; i++) {
        rays.push_back(returnGlossyRay(reflection, degreeBlur));
    }
    return rays;
}

glm::vec3 getPlaneCoord(Ray ray, glm::vec3 point) {
    glm::vec3 w = glm::normalize(ray.direction);
    glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
    glm::vec3 x = glm::normalize(glm::cross(t, w));
    glm::vec3 y = glm::normalize(glm::cross(w, x));
    glm::vec3 rToP = point - ray.origin;
    return ray.origin + glm::dot(rToP, x) * x + glm::dot(rToP, y) * y;
}

glm::vec3 trilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features, const Ray& ray, HitInfo hitInfo)
{
    ImageMipMap mipmap = getMipMap(image);
    glm::vec3 point = ray.origin + ray.t * ray.direction;
    float dist = glm::distance(point, ray.origin);
    float angle = acos(glm::dot(-ray.direction, hitInfo.normal));

    Mesh& mesh = *hitInfo.mesh;
    glm::uvec3 tri = hitInfo.triangle;
    // If a sphere was hit first, we compute bilinear interpolation
    if (tri == glm::uvec3(-1000000)) {
        return bilinearInterpolation(image, texCoord, features);
    }
    glm::vec3 v0 = mesh.vertices[tri[0]].position;
    glm::vec3 v1 = mesh.vertices[tri[1]].position;
    glm::vec3 v2 = mesh.vertices[tri[2]].position;

    glm::vec3 p0 = getPlaneCoord(ray, v0);
    glm::vec3 p1 = getPlaneCoord(ray, v1);
    glm::vec3 p2 = getPlaneCoord(ray, v2);

    float areaV = glm::length(glm::cross(v0 - v1, v0 - v2));
    float areaP = glm::length(glm::cross(p0 - p1, p0 - p2));

    float k = (areaV / areaP - 1.0f) / 2.8f;

    //std::cout << k << " " << areaV << " " << areaP << std::endl;

    //float k = dist * angle / 3.0f;
    glm::vec3 w = glm::normalize(hitInfo.normal);
    glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
    glm::vec3 xVector = glm::normalize(glm::cross(t, w));
    glm::vec3 yVector = glm::normalize(glm::cross(w, xVector));
    Ray rayU = Ray { point, yVector, areaV / 5.0f };
    Ray rayV = Ray { point, xVector, areaP / 5.0f };
    drawRay(rayU, glm::vec3 { 1.0f, 0.5f, 0.0f });
    drawRay(rayV, glm::vec3 { 1.0f, 0.0f, 0.5f });
    float k0 = std::floor(k);
    float k1 = k0 + 1;
    float a = k1 - k;
    if (showMipmapLevel == true) {
        int thisK = mipmapLevel;
        if (thisK > mipmap.height.size() - 1) {
            thisK = mipmap.height.size() - 1;
        }
        return bilinearInterpolationForMipMap(mipmap, thisK, texCoord, features);
    }
    if (k0 < 0) {
        return bilinearInterpolation(image, texCoord, features);
    }
    if (k1 >= mipmap.height.size()) {
        return bilinearInterpolationForMipMap(mipmap, mipmap.height.size() - 1, texCoord, features);
    }
    glm::vec3 c0 = bilinearInterpolationForMipMap(mipmap, k0, texCoord, features);
    glm::vec3 c1 = bilinearInterpolationForMipMap(mipmap, k1, texCoord, features);
    return a * c0 + (1 - a) * c1;
}

glm::vec3 bilinearInterpolationForMipMap(const ImageMipMap& image, int level, const glm::vec2& texCoord, const Features& features)
{
    glm::vec2 texelPos { (image.width[level] - 1) * texCoord[0], (image.height[level] - 1) * (1 - texCoord[1]) };
    //texelPos.x = std::max(0.0f, std::min(float(image.width[level] - 1), texelPos.x));
    //texelPos.y = std::max(0.0f, std::min(float(image.height[level] - 1), texelPos.y));
    glm::vec2 lowerPos { floor((image.width[level] - 1) * texCoord[0]), floor((image.height[level] - 1) * (1 - texCoord[1])) };
    glm::vec2 upperPos { lowerPos.x + 1, lowerPos.y + 1 };
    float u = texelPos.x - lowerPos.x;
    float v = texelPos.y - lowerPos.y;
    lowerPos = glm::mod(lowerPos, glm::vec2 { image.width[level], image.height[level] });
    upperPos = glm::mod(upperPos, glm::vec2 { image.width[level], image.height[level] });
    glm::vec3 lowerLeft = image.pixels[level][lowerPos.y * image.width[level] + lowerPos.x];
    glm::vec3 upperRight = image.pixels[level][upperPos.y * image.width[level] + upperPos.x];
    glm::vec3 lowerRight = image.pixels[level][upperPos.y * image.width[level] + lowerPos.x];
    glm::vec3 upperLeft = image.pixels[level][lowerPos.y * image.width[level] + upperPos.x];
    return lowerLeft * (1.0f - u) * (1.0f - v) + upperRight * u * v + lowerRight * (1.0f - u) * v + upperLeft * u * (1.0f - v);
}