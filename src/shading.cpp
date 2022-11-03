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
    glm::vec3 view = ray.origin - position;
    glm::vec3 h = glm::normalize(lightDir + view);
    if (glm::dot(h, glm::normalize(hitInfo.normal)) < 0.0f) {
        return {0, 0, 0};
    }
    // EXPERIMENT
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
        return lightColor * texel * dot;
    }
    return lightColor * hitInfo.material.kd * dot;
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

glm::vec3 trilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features, const Ray& ray, HitInfo hitInfo)
{
    // Get the mipmap first
    
    ImageMipMap mipmap = getMipMap(image);
    
    /** for (int i = 0; i < mipmap.height.size(); i++) {
        std::cout << "width: " << mipmap.width[i] << " height: " << mipmap.height[i] << " size: " << mipmap.pixels[i].size() << std::endl;
    }**/
    glm::vec3 point = ray.origin + ray.t * ray.direction;
    Ray rayCopy = ray;
    // Get the triangle which was hit, in order to calculate the texCoord and 3D coordinates for a corner
    //Mesh& mesh = *hitInfo.mesh;
    //glm::uvec3 tri = hitInfo.triangle;
    // If a sphere was hit first, we compute bilinear interpolation
    //if (tri == glm::uvec3(-1000000)) {
      //  return bilinearInterpolation(image, texCoord, features);
    //}
    //auto v = mesh.vertices[tri[0]];
    //if (point == v.position) {
      //  v = mesh.vertices[tri[1]];
    //}
    // Computing the derivative
    // First we calculate the 2D distance between the texture coordinates
    /** glm::vec2 texCoord2 = v.texCoord;
    float distanceTexX = fabs(texCoord.x - texCoord2.x);
    float distanceTexY = fabs(texCoord.y - texCoord2.y);
    // Then we approximate the screen space
    Ray newRay = ray;
    newRay.direction = glm::normalize(v.position - ray.origin) * glm::length(ray.direction);
    newRay.t = ray.t;
    glm::vec3 newPoint = newRay.t * newRay.direction + newRay.origin;
    float distanceScreenX = fabs(point.x - newPoint.x);
    float distanceScreenY = fabs(point.y - newPoint.y);
    // The derivative is the distanceScreen / distanceTex
    float derivativeX = 1.0, derivativeY = 1.0;
    // If the distance for coordinates is zero, the derivative will converge to 1, so it stays 1
    if (distanceScreenX != 0.0f) {
        derivativeX = distanceTexX / distanceScreenX;
    }
    // Same for y axis
    if (distanceScreenY != 0.0f) {
        derivativeY = distanceTexY / distanceScreenY;
    }
    // Taking the maximum derivative
    float maxDerivative = derivativeX;
    if (maxDerivative > derivativeY) {
        maxDerivative = derivativeY;
    }**/
    float dist = glm::distance(point, ray.origin);
    float angle = acos(glm::dot(-ray.direction, hitInfo.normal));
        // Finding the log
        // float k = std::log2(maxDerivative);
    float k = dist * angle / 3.0f;
    glm::vec3 w = glm::normalize(hitInfo.normal);
    glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
    glm::vec3 xVector = glm::normalize(glm::cross(t, w));
    glm::vec3 yVector = glm::normalize(glm::cross(w, xVector));
    Ray rayU = Ray { point, yVector, dist / 5.0f };
    Ray rayV = Ray { point, xVector, angle / 5.0f };
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
    //std::cout << k0 << " " << k << " " << k1 << " Distance: " << dist << " Angle: " << angle << std::endl;
    // The level cannot be negative, so we approximate the trilinear interpolation to a bilinear interpolation for the first level in the mipmap (level 0)
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
    texelPos.x = std::max(0.0f, std::min(float(image.width[level] - 1), texelPos.x));
    texelPos.y = std::max(0.0f, std::min(float(image.height[level] - 1), texelPos.y));
    glm::vec2 lowerPos { floor((image.width[level] - 1) * texCoord[0]), floor((image.height[level] - 1) * (1 - texCoord[1])) };
    glm::vec2 upperPos { lowerPos.x + 1, lowerPos.y + 1 };
    float u = texelPos.x - lowerPos.x;
    float v = texelPos.y - lowerPos.y;
    //std::cout << texelPos.x << " " << texelPos.y << std::endl;
    lowerPos = glm::mod(lowerPos, glm::vec2 { image.width[level], image.height[level] });
    upperPos = glm::mod(upperPos, glm::vec2 { image.width[level], image.height[level] });
    glm::vec3 upperLeft = image.pixels[level][lowerPos.y * image.width[level] + lowerPos.x];
    glm::vec3 lowerRight = image.pixels[level][upperPos.y * image.width[level] + upperPos.x];
    glm::vec3 upperRight = image.pixels[level][upperPos.y * image.width[level] + lowerPos.x];
    glm::vec3 lowerLeft = image.pixels[level][lowerPos.y * image.width[level] + upperPos.x];
    return upperLeft * (1.0f - u) * (1.0f - v) + lowerRight * u * v + upperRight * (1.0f - u) * v + lowerLeft * u * (1.0f - v);
}