#include "texture.h"
#include <framework/image.h>
#include <iostream>
#include <glm/gtx/common.hpp>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    int i = (image.width * texCoord[0]);
    int j = (image.height * texCoord[1]);

    // clamp if i and j are out of range

    i = std::max(0, std::min(image.width - 1, i));
    j = std::max(0, std::min(image.height - 1, j));
    return image.pixels[j*image.width + i];
}

glm::vec3 bilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features) {
    glm::vec2 texelPos { (image.width - 1) * texCoord[0], (image.height - 1) * (1 - texCoord[1])};
    texelPos.x = std::max(0.0f, std::min(float(image.width - 1), texelPos.x));
    texelPos.y = std::max(0.0f, std::min(float(image.height - 1), texelPos.y));
    glm::vec2 lowerPos { floor((image.width - 1) * texCoord[0]), floor((image.height - 1) * texCoord[1]) };
    glm::vec2 upperPos { lowerPos.x + 1, lowerPos.y + 1 };
    float u = texelPos.x - lowerPos.x;
    float v = texelPos.y - lowerPos.y;
    //lowerPos = glm::mod(lowerPos, glm::vec2 { image.width, image.height });
    //upperPos = glm::mod(upperPos, glm::vec2 { image.width, image.height });
    glm::vec3 upperLeft = image.pixels[lowerPos.y * image.width + lowerPos.x];
    glm::vec3 lowerRight = image.pixels[upperPos.y * image.width + upperPos.x];
    glm::vec3 upperRight = image.pixels[upperPos.y * image.width + lowerPos.x];
    glm::vec3 lowerLeft = image.pixels[lowerPos.y * image.width + upperPos.x];
    return upperLeft * (1.0f - u) * (1.0f - v) + lowerRight * u * v + upperRight * (1.0f - u) * v + lowerLeft * u * (1.0f - v);
}

glm::vec2 getUVForBilinear(const Image& image, const glm::vec2& texCoord, const Features& features) {
    glm::vec2 texelPos { (image.width - 1) * texCoord[0], (image.height - 1) * texCoord[1] };
    texelPos.x = std::max(0.0f, std::min(float(image.width - 1), texelPos.x));
    texelPos.y = std::max(0.0f, std::min(float(image.height - 1), texelPos.y));
    glm::vec2 lowerPos { floor((image.width - 1) * texCoord[0]), floor((image.height - 1) * texCoord[1]) };
    glm::vec2 upperPos { lowerPos.x + 1, lowerPos.y + 1 };
    float u = texelPos.x - lowerPos.x;
    float v = texelPos.y - lowerPos.y;
    return glm::vec2 { u, v };
}
