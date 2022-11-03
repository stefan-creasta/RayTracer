#include "mrays.h"
#include "draw.h"
#include <random>
#include "render.h"



glm::vec3 calculateColorMultipleRaysPerPixel(const Scene &scene, Ray &ray, const BvhInterface& bvh,  const Features &features, int rayDepth, int samples) {

    std::default_random_engine rng;
    std::uniform_real_distribution<float> dist(-0.1, 0.1);

    glm::vec3 lower = ray.origin - glm::vec3{0.1,0.1,0};
    glm::vec3 upper = ray.origin + glm::vec3{0.1,0.1,0};
    drawAABB(AxisAlignedBox({lower, upper}));

    glm::vec3 fcolor(0.0f);
    for (int i = 0; i < samples; i++) {
        glm::vec3 r = ray.origin + glm::vec3{dist(rng), dist(rng), 0};
        Ray temp = {r, ray.direction};
        auto col = getFinalColor(scene, bvh, temp, features, 0);
        fcolor += col;

    }

    fcolor *= glm::vec3(1.0f/(samples*1.0f)); 

   /* for (int i = 0; i < samples; i++) {
        float xeps = dist(rng) / windowResolution.x * 2.0f;
        float yeps = dist(rng) / windowResolution.y * 2.0f;

        Ray r = {ray.origin + {xeps, }}
    }*/

    return fcolor;
}