#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <random>
bool flag = false; // debug random light samples
#include <iostream>
glm::vec3 posForDebug = glm::vec3(-10000.0f);
std::vector<glm::vec3> prevPos;
std::vector<glm::vec3> prevCol;


float getRandomVal()
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float randomVal = distribution(generator);
    return randomVal;
}
// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    glm::vec3 p0ToP1 = segmentLight.endpoint1 - segmentLight.endpoint0;
    float alpha = getRandomVal();
    position = alpha * p0ToP1 + segmentLight.endpoint0;
    color = (1 - alpha) * segmentLight.color0 + alpha * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    // first v0 to v1
    float alpha1 = getRandomVal();
    glm::vec3 c1 = (1 - alpha1) * parallelogramLight.color0 + alpha1 * parallelogramLight.color1;
    glm::vec3 c2 = (1 - alpha1) * parallelogramLight.color2 + alpha1 * parallelogramLight.color3;
    float alpha2 = getRandomVal();
    color = (1 - alpha2) * c1 + alpha2 * c2;
    position = parallelogramLight.v0 + parallelogramLight.edge01 * alpha1 + parallelogramLight.edge02 * alpha2;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 hit = ray.origin + ray.t * ray.direction;
    float eps = 1e-6;
    glm::vec3 normal = hitInfo.normal;
    if (glm::dot(glm::normalize(hitInfo.normal), glm::normalize(ray.origin - hit)) < -eps) {
        normal = -hitInfo.normal;
    }
    
    float length = glm::length(hit - samplePos);
    glm::vec3 copyPos = samplePos;
    Ray sray { copyPos, hit - samplePos };
    HitInfo newhit;
    bvh.intersect(sray, newhit, features);
    
    glm::vec3 secondHit = sray.origin + sray.t * sray.direction;
    if (glm::dot(glm::normalize(samplePos - hit), glm::normalize(normal)) < -eps) {

        if ((features.enableHardShadow || features.enableSoftShadow))
            drawRay(sray, glm::vec3{1, 0, 0});
        return 0.0f;
    }
    if (glm::distance(hit, secondHit) > 1e-3) {
        if ((features.enableHardShadow || features.enableSoftShadow))
            drawRay(sray, glm::vec3{1, 0, 0});
        return 0.0f;
    }

    if ((features.enableHardShadow || features.enableSoftShadow))
        drawRay(sray, debugColor);
    
    return 1.0f;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.


glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        glm::vec3 med = {0, 0, 0};

        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                if (features.enableHardShadow)
                    med += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo) * testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo);
                else
                    med += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);

            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                // Perform your calculations for a segment light.
                if (features.enableSoftShadow) {
                    int sampleSize = 50;
                    glm::vec3 avgColor = { 0.0f, 0.0f, 0.0f };
                    //DEBUG FOR SEGMENT LIGHT
                    glm::vec3 pos = ray.t * ray.direction + ray.origin;
                    float eps = 1e-6;
                    if (std::fabs(posForDebug.x - pos.x) < eps || std::fabs(posForDebug.y - pos.y) < eps || std::fabs(posForDebug.z - pos.z) < eps) {
                        flag = false;
                        prevCol.clear();
                        prevPos.clear();
                    }
                    //DEBUG FOR SEGMENT LIGHT
                    for (int i = 0; i < sampleSize; i++) {
                        glm::vec3 pos;
                        glm::vec3 col;
                        if (flag == false) {
                            sampleSegmentLight(segmentLight, pos, col);
                            prevPos.push_back(pos);
                            prevCol.push_back(col);
                        } else {
                            pos = prevPos[i];
                            col = prevCol[i];
                        }
                        avgColor += computeShading(pos, col, features, ray, hitInfo) * testVisibilityLightSample(pos, col, bvh, features, ray, hitInfo);
                    }
                    flag = true;
                    med += avgColor * float((1.0 / float(sampleSize)));
                }
            } else if (std::holds_alternative<ParallelogramLight>(light)) { 
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                // Perform your calculations for a parallelogram light.
                if (features.enableSoftShadow) {
                    int sampleSize = 50;
                    glm::vec3 avgColor = { 0.0f, 0.0f, 0.0f };
                    // DEBUG FOR PARALLELOGRAM LIGHT
                    glm::vec3 pos = ray.t * ray.direction + ray.origin;
                    float eps = 1e-6;
                    if (std::fabs(posForDebug.x - pos.x) < eps || std::fabs(posForDebug.y - pos.y) < eps || std::fabs(posForDebug.z - pos.z) < eps) {
                        flag = false;
                        prevCol.clear();
                        prevPos.clear();
                    }
                    // DEBUG FOR PARALLELOGRAM LIGHT
                    for (int i = 0; i < sampleSize; i++) {
                        glm::vec3 pos;
                        glm::vec3 col;
                        if (flag == false) {
                            sampleParallelogramLight(parallelogramLight, pos, col);
                            prevPos.push_back(pos);
                            prevCol.push_back(col);
                        } else {
                            pos = prevPos[i];
                            col = prevCol[i];
                        }
                        // std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;
                        avgColor += computeShading(pos, col, features, ray, hitInfo) * testVisibilityLightSample(pos, col, bvh, features, ray, hitInfo);
                    }
                    flag = true;
                    med += avgColor * float((1.0 / float(sampleSize)));
                }
            }
        }
        return med;
    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
