#include <framework/image.h>
#include <framework/ray.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
#include <glm/trigonometric.hpp>
#include <glm/vector_relational.hpp>
#include <glm/common.hpp>
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <optional>
#include "environment_mapping.h"
#include "texture.h"
#include <iostream>

EnvironmentMap EnvironmentMap::loadEnvironmentMap(const std::filesystem::path& file, EnvironmentMappingType mappingType, float verticalFieldOfView, glm::vec3 backgroundColor)
{
    Image image(file);
    return EnvironmentMap(image, mappingType, verticalFieldOfView, backgroundColor);
}

EnvironmentMap::EnvironmentMap(glm::vec3 color)
    : backgroundColor(color)
{
}

EnvironmentMap::EnvironmentMap(Image& texture, EnvironmentMappingType mappingType, float verticalFieldOfView, glm::vec3 backgroundColor)
    : backgroundColor(backgroundColor)
    , texture(texture)
    , verticalFOVFactor(glm::radians(verticalFieldOfView))
    , mappingType(mappingType)
{
}
#define PI 3.14159265358979323846
glm::vec3 EnvironmentMap::getColor(Ray ray, const Features& features) const
{
    if (texture.has_value()) {
        float x, y;

        x = glm::atan(ray.direction.z, ray.direction.x) / (2 * PI) + 0.5;
        switch (this->mappingType) {
        case SPHERICAL:
            y = (PI / this->verticalFOVFactor) * (glm::acos(glm::dot(glm::normalize(ray.direction), glm::vec3 { 0.f, 1.f, 0.f })) / PI - 0.5) + 0.5;
            break;
        case CYLINDRICAL:
            y = -0.5 / glm::tan(0.5f * this->verticalFOVFactor) / glm::tan(glm::acos(glm::dot(glm::normalize(ray.direction), glm::vec3 { 0.f, 1.f, 0.f }))) + 0.5;
            break;
        case SPHEROCYLINDRICAL:
            y = -0.5 * PI / this->verticalFOVFactor * glm::dot(glm::normalize(ray.direction), glm::vec3 { 0.f, 1.f, 0.f }) + 0.5;
            break;
        default:
            y = 0.0;
        }
        //const float cosy = glm::dot(glm::normalize(ray.direction), glm::vec3 { 0.f, 1.f, 0.f });
        //const float unadjustedY = glm::acos(cosy) / PI;
        //const float y = 0.5f * (this->verticalFOVFactor * cosy + 1.0);
        //const float y = this->verticalFOVFactor * ()
        if (glm::abs(y - 0.5) > 0.5)
            return this->backgroundColor;

        if (features.extra.enableBilinearTextureFiltering) { 
            return bilinearInterpolation(texture.value(), glm::vec2 {x, y}, features);
        } else {
            return acquireTexel(texture.value(), glm::vec2 { x, y }, features);
        }
    } else {
        return this->backgroundColor;
    }
}
