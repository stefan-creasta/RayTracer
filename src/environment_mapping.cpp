#pragma once
#include <framework/ray.h>
DISABLE_WARNINGS_PUSH()
#include <glm/trigonometric.hpp>
#include <glm/vector_relational.hpp>
#include <glm/common.hpp>
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include "environment_mapping.h"
#include "texture.h"
#include <iostream>
#include <random>
#include "sampling.h"

#define PI 3.14159265358979323846

float getRandomValForEnvironmentMapping()
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float randomVal = distribution(generator);
    return randomVal;
}

float getRadiance(glm::vec3 color)
{
    return 0.25 * (color.x + 2. * color.y + color.z);
}

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
    float totalRadiance = 0.f;
    for (const glm::vec3& pixel : texture.pixels) {
        totalRadiance += getRadiance(pixel);
    }
    buildRadianceHierarchy({}, totalRadiance);
}

void EnvironmentMap::buildRadianceHierarchy(const AxisAlignedRectangle& aar, const float totalRadiance)
{
    const glm::vec2 aarSize = aar.upper - aar.lower;

    const float minRadiance = (texture.value().width * texture.value().height) / 1000;
    if (!texture || totalRadiance < minRadiance || (aarSize.x + aarSize.y) < 0.1) {
        radianceBins.push_back(aar);
        return;
    }

    const glm::vec2 imageSize = glm::vec2 { texture.value().width, texture.value().height };
    const glm::vec2 imageSizeM1 = glm::vec2 { texture.value().width - 1, texture.value().height - 1 };

    bool splitVertically = aarSize.y / aarSize.x > 1.f;
    //std::cout << splitVertically << std::endl;
    const glm::vec2 splitUpper = splitVertically
        ? glm::vec2 {aar.upper.x, 0.5 * (aar.lower.y + aar.upper.y)}
        : glm::vec2 {0.5 * (aar.lower.x + aar.upper.x), aar.upper.y};

    const glm::vec2 lowerScaled = aar.lower * imageSizeM1;
    const glm::vec2 splitUpperScaled = splitUpper * imageSizeM1;

    float lowerHalfRadiance = 0.f;
    for (int i = lowerScaled.x; i < splitUpperScaled.x; i++) {
        for (int j = lowerScaled.y; j < splitUpperScaled.y; j++) { 
            lowerHalfRadiance += getRadiance(texture.value().pixels[i + j * imageSize.x]);
        }
    }

    float upperHalfRadiance = totalRadiance - lowerHalfRadiance;

    const AxisAlignedRectangle aarLower 
    {
        aar.lower,
        splitUpper
    };

    const AxisAlignedRectangle aarUpper {
        splitVertically ? glm::vec2 { aar.lower.x, splitUpper.y } : glm::vec2 { splitUpper.x, aar.lower.y },
        aar.upper
    };

    if (lowerHalfRadiance > 0.5 * totalRadiance)
        buildRadianceHierarchy(aarLower, lowerHalfRadiance);
    else
        radianceBins.push_back(aarLower);
    if (upperHalfRadiance > 0.5 * totalRadiance)
        buildRadianceHierarchy(aarUpper, upperHalfRadiance);
    else
        radianceBins.push_back(aarUpper);
}

std::vector<Ray> EnvironmentMap::getSamplingRay(const glm::vec3& position, const glm::vec3& normal, int n) const
{
    if (texture) {
        int nPerTry = glm::max(glm::sqrt(n / 10.), 1.);
        std::vector<Ray> returned;
        returned.reserve(n);
        for (int i = 0; returned.size() < n && i < 100; i++) {
            const int randomBin = std::floor(getRandomValForEnvironmentMapping() * radianceBins.size());
            const AxisAlignedRectangle& rect = radianceBins[randomBin];
            std::vector<glm::vec2> samples = sample2D(rect, nPerTry, nPerTry);
            for (glm::vec2& sample : samples) { 
                Ray ray = getRayForCoordinate(sample);
                const float eps = 0.0001f / glm::dot(glm::normalize(ray.direction), normal);
                ray.origin = position + eps * ray.direction;
                if (glm::dot(ray.direction, normal) > 0.f)
                    returned.push_back(ray);
            }
        }
        return returned;
    }
    return { Ray { position, normal } };
}

Ray EnvironmentMap::getRayForCoordinate(const glm::vec2& coordinates) const { 
    const float azimuth = (2.f * PI * (coordinates.x - 0.5));
    float pitch;
    glm::vec3 pitched, yawed;
    switch (mappingType) {
    case SPHERICAL:
        pitch = (0.5 - coordinates.y) * verticalFOVFactor;
        pitched = {
            glm::cos(pitch),
            glm::sin(pitch),
            0.f,
        };
        yawed = {
            pitched.x * glm::cos(azimuth),
            pitched.y,
            pitched.x * glm::sin(azimuth),
        };
        return Ray { {}, yawed };
        break;
    case CYLINDRICAL:
        pitch = glm::atan(glm::tan(0.5f * this->verticalFOVFactor) * (0.5 - coordinates.y) / 0.5);
        pitched = {
            glm::cos(pitch),
            glm::sin(pitch),
            0.f,
        };
        yawed = {
            pitched.x * glm::cos(azimuth),
            pitched.y,
            pitched.x * glm::sin(azimuth),
        };
        return Ray { {}, yawed };
        break;
    }
    return Ray { {}, { 1.f, 0.f, 0.f } };
}

glm::vec3 EnvironmentMap::getColor(Ray ray, const Features& features) const
{
    if (texture.has_value()) {
        float x, y;

        x = glm::atan(ray.direction.z, ray.direction.x) / (2 * PI) + 0.5;
        switch (this->mappingType) {
        case SPHERICAL:
            y = -(PI / this->verticalFOVFactor) * (glm::acos(glm::dot(glm::normalize(ray.direction), glm::vec3 { 0.f, 1.f, 0.f })) / PI - 0.5) + 0.5;
            break;
        case CYLINDRICAL:
            y = 0.5 / glm::tan(0.5f * this->verticalFOVFactor) / glm::tan(glm::acos(glm::dot(glm::normalize(ray.direction), glm::vec3 { 0.f, 1.f, 0.f }))) + 0.5;
            break;
        }

        if (glm::abs(y - 0.5) > 0.5)
            return this->backgroundColor;

        // Debug: Show the Radiance Bins
        /*for (const AxisAlignedRectangle& rectangle : radianceBins) { 
            glm::vec2 middle = 0.5f * (rectangle.lower + rectangle.upper);
            middle.y = 1.f - middle.y;
            const glm::vec2 aabSize = rectangle.upper - rectangle.lower;
            const float diffx = 0.5 * aabSize.x - glm::abs(x - middle.x);
            const float diffy = 0.5 * aabSize.y - glm::abs(y - middle.y);

            if (diffx > 0 && diffy > 0 && (diffx < 0.01 || diffy < 0.01))
                return glm::vec3 { 1.f, 0.f, 0.f };
        }*/

        if (features.extra.enableBilinearTextureFiltering) { 
            return bilinearInterpolation(texture.value(), glm::vec2 {x, y}, features);
        } else {
            return acquireTexel(texture.value(), glm::vec2 { x, y }, features);
        }
    } else {
        return this->backgroundColor;
    }
}
