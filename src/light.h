#pragma once
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "shading.h"

extern int sampleSize;

void sampleSegmentLight (const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color);

void sampleParallelogramLight (const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color);

glm::vec3 sampleEnvironment(const EnvironmentMap& map, const BvhInterface& bvh, const Ray& ray, const HitInfo& hitInfo, const Features& features);

float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

