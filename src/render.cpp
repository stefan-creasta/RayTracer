#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "transparency.h"
#include "dof.h"
#include "environment_mapping.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>
#include <random>


glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {
        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
        glm::vec3 color = { 1, 1, 1 };
        if (features.enableRecursive) { 
            std::vector<Ray> rays;
            Ray reflection = computeReflectionRay(ray, hitInfo);
            if (!features.extra.enableGlossyReflection) {
                rays.push_back(reflection);
            } else {
                rays = glossyRays(reflection);
            }
            if (rayDepth > 0 && hitInfo.material.ks != glm::vec3 {0.0, 0.0, 0.0}) {
                color = { 0, 0, 0 };
                for (int i = 0; i < rays.size(); i++) {
                    color += getFinalColor(scene, bvh, rays[i], features, rayDepth - 1);
                }
                color *= (1.0f / float(rays.size()));
                
            }
            else {
                drawRay(ray, Lo);
                return Lo;
            }
            drawRay(ray, color * hitInfo.material.ks + Lo);
            return (Lo + color * hitInfo.material.ks);
        }
        drawRay(ray, Lo);
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        
        // Set the color of the pixel to the environment color (black if environment mapping is disabled) if the ray misses.
        if (features.extra.enableEnvironmentMapping) {
            glm::vec3 color = (* scene.environmentMap[0]).getColor(ray, features);
            drawRay(ray, color);
            return color;
        } else {
            drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
            return glm::vec3 { 0.f };
        }
    }
}


void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode

    float aperture = 0.1;
    std::default_random_engine rng;
    std::uniform_real_distribution<float> dist(-aperture/2, aperture/2);
    #ifdef NDEBUG
    #pragma omp parallel for schedule(guided)
    #endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            Ray cameraRay = camera.generateRay(normalizedPixelPos);
            /*
                Extra Feature: Transparency
            */
            glm::vec3 col;
            col = getFinalColor(scene, bvh, cameraRay, features, 1);
            screen.setPixel(x, y, col);
            
        }
    }
}

void renderRayTracingDepthOfField(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features, float aperture, float focalLength, int samples)
{
    glm::ivec2 windowResolution = screen.resolution();

    std::default_random_engine rng;
    std::uniform_real_distribution<float> dist(-aperture/2, aperture/2);

    // Enable multi threading in Release mode
    #ifdef NDEBUG
    #pragma omp parallel for schedule(guided)
    #endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            Ray cameraRay = camera.generateRay(normalizedPixelPos);
            glm::vec3 col = glm::vec3 {0, 0, 0};
                
            glm::vec3 focalPoint = getFocalPoint(cameraRay, focalLength);
            std::vector <Ray> rays = sampledRays(focalPoint, generateSamples(cameraRay.origin, aperture, samples, rng, dist));
            for (auto r : rays) {
                col += getFinalColor(scene, bvh, r, features, 0);
            }
            col = col * glm::vec3(1.0f / (rays.size()));
            screen.setPixel(x, y, col);
        }
    }
}


void renderRayTracingTransparency(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
    #ifdef NDEBUG
    #pragma omp parallel for schedule(guided)
    #endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            Ray cameraRay = camera.generateRay(normalizedPixelPos);
            glm::vec3 col;
            col = calculateColorTransparency(scene, cameraRay, bvh, features, 1);
            screen.setPixel(x, y, col);
            
        }
    }
}

void renderRayTracingMRaysPerPixel(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features, int samples) {
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
    

    std::default_random_engine rng;
    std::uniform_real_distribution<float> dist(0, 1);
    #ifdef NDEBUG
    #pragma omp parallel for schedule(guided)
    #endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            
            glm::vec3 fcolor(0.0f);
            for (int i = 0; i < samples; i++) {
                const glm::vec2 normalizedPixelPos {
                    float(x + dist(rng)) / float(windowResolution.x) * 2.0f - 1.0f,
                    float(y + dist(rng)) / float(windowResolution.y) * 2.0f - 1.0f
                };
                Ray cameraRay = camera.generateRay(normalizedPixelPos);
                fcolor += getFinalColor(scene, bvh, cameraRay, features, 0);
            }
            fcolor *= glm::vec3(1.0f/(samples*1.0f));
            screen.setPixel(x, y, fcolor);
        }
    }
}

glm::vec3 startLookAt;
glm::vec3 startAngles;
float startDistance;
bool rendered;

void renderRayTracingMotionBlur(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features, int steps)
{
    glm::ivec2 windowResolution = screen.resolution();
    std::default_random_engine rng;
    std::uniform_real_distribution<float> rand(0.f, 1.f);

    if (!rendered) { 
        renderRayTracing(scene, camera, bvh, screen, features);
        startLookAt = camera.lookAt();
        startAngles = camera.rotationEulerAngles();
        startDistance = camera.distanceFromLookAt();
        rendered = true;
        return;
    }

    glm::vec3 endLookAt = camera.lookAt();
    glm::vec3 endAngles = camera.rotationEulerAngles();
    float endDistance = camera.distanceFromLookAt();

    glm::vec3 scaleLookAt = endLookAt - startLookAt;
    glm::vec3 scaleAngles = endAngles - startAngles;
    float scaleDistance = endDistance - startDistance;

    startLookAt = camera.lookAt();
    startAngles = camera.rotationEulerAngles();
    startDistance = camera.distanceFromLookAt();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            glm::vec3 col {0.f};
            for (int i = 0; i < steps; i++) {
                float jitter = rand(rng);
                float lerpFactor = jitter + i;
                glm::vec3 lerpLookAt = startLookAt + (lerpFactor / (steps)) * scaleLookAt;
                glm::vec3 lerpAngles = startAngles + (lerpFactor / (steps)) * scaleAngles;
                float lerpDistance = startDistance + (lerpFactor / (steps)) * scaleDistance;

                Trackball camera2 = camera;
                camera2.setCamera(lerpLookAt, lerpAngles, lerpDistance);

                // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
                const glm::vec2 normalizedPixelPos {
                    float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                    float(y) / float(windowResolution.y) * 2.0f - 1.0f
                };
                Ray cameraRay = camera2.generateRay(normalizedPixelPos);
                /*
                    Extra Feature: Transparency
                */
                
                col += getFinalColor(scene, bvh, cameraRay, features, 1);
            }
            screen.setPixel(x, y, col / (float) steps);
        }
    }
}