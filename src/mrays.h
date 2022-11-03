#pragma once


#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include "common.h"
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()

#include "bvh_interface.h"

struct Ray;
struct Scene;
struct Features;

/*

    EXTRA FEATURE: MULTIPLE RAYS PER PIXEL WITH IRREGULAR SAMPLING
    IMPLEMENTATION: mrays.cpp

*/

glm::vec3 calculateColorMultipleRaysPerPixel(const Scene &scene, Ray &ray, const BvhInterface& bvh,  const Features &features, int rayDepth, int samples = 64);