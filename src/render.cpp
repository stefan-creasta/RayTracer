#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "transparency.h"
#include "dof.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>
#include <random>

std::vector <Ray> defaultRaysDOF;
int rD = -1;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {
        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
        glm::vec3 color = { 1, 1, 1 };
        if (features.enableRecursive) { // HARD STOP DEPTH OF FIELD
            Ray reflection = computeReflectionRay(ray, hitInfo);
            if (rayDepth > 0 && hitInfo.material.ks != glm::vec3 {0.0, 0.0, 0.0}) {
                color = getFinalColor(scene, bvh, reflection, features, rayDepth - 1);
            }
            else {
                drawRay(ray, Lo);
                return Lo * hitInfo.material.transparency;
            }
            drawRay(ray, color * hitInfo.material.ks + Lo);
            return (Lo + color * hitInfo.material.ks) * hitInfo.material.transparency;
        }
        drawRay(ray, Lo);
        return Lo * hitInfo.material.transparency;
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

    float aperture = 0.5;
    std::default_random_engine rng{42};
    std::uniform_real_distribution<float> dist(-aperture/2, aperture/2);
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
            if (features.extra.enableTransparency) {
                col = calculateColorTransparency(scene, cameraRay, bvh, features, 1);
            }
            else {
                col = getFinalColor(scene, bvh, cameraRay, features, 1);
            }
            if (features.extra.enableDepthOfField) {

                
                col = glm::vec3 {0, 0, 0};
                float focalLength = 1.5;
                int samples = 100;
                
                glm::vec3 focalPoint = getFocalPoint(cameraRay, focalLength);
                std::vector <Ray> rays = sampledRays(focalPoint, generateSamples(cameraRay.origin, aperture, samples, rng, dist));
                for (auto r : rays) {
                    col += getFinalColor(scene, bvh, r, features, 0);
                    std::cout << col[0] << " " << col[1] << " " << col[2] << "\n"; 
                }

                col = col * glm::vec3(1 / (rays.size()));
                //std::cout << col[0] << " " << col[1] << " " << col[2] << "\n"; 
            }
            screen.setPixel(x, y, col);
            
        }
    }
}