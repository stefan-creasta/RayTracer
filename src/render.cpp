#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {
        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
        glm::vec3 color = { 1, 1, 1 };
        if (features.enableRecursive) { 
            std::vector<Ray> rays;
            if (!features.extra.enableGlossyReflection) {
                rays.push_back(computeReflectionRay(ray, hitInfo));
            } else {
                rays = glossyRays(ray, hitInfo);
            }
            if (rayDepth > 0 && hitInfo.material.ks != glm::vec3 {0.0, 0.0, 0.0}) {
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
            return Lo + color * hitInfo.material.ks;
        }
        drawRay(ray, Lo);
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
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
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features, 1));
        }
    }
}