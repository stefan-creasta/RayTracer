#include "texture.h"
#include "draw.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <random>
#include <draw.cpp>
#include <iostream>
float degreeBlur = 0.01f;
int numberOfRays = 5;

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

Ray returnGlossyRay(Ray reflection)
{
    glm::vec3 w = glm::normalize(reflection.direction);
    glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
    glm::vec3 u = glm::normalize(glm::cross(t, w));
    glm::vec3 v = glm::normalize(glm::cross(w, u));
    float ua = -degreeBlur / 2.0f + degreeBlur * getRandomVal2();
    float va = -degreeBlur / 2.0f + degreeBlur * getRandomVal2();
    Ray rr = reflection;
    rr.direction += ua * u + va * v;
    glm::vec3 p1 = reflection.origin + u * degreeBlur + v * degreeBlur;
    glm::vec3 p2 = reflection.origin + u * degreeBlur - v * degreeBlur;
    glm::vec3 p3 = reflection.origin - u * degreeBlur + v * degreeBlur;
    glm::vec3 p4 = reflection.origin - u * degreeBlur - v * degreeBlur;
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

std::vector<Ray> glossyRays(Ray reflection) {
    std::vector<Ray> rays;
    for (int i = 1; i <= numberOfRays; i++) {
        rays.push_back(returnGlossyRay(reflection));
    }
    return rays;
}