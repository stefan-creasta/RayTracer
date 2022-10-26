#include "texture.h"
#include <framework/image.h>
#include <iostream>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    int i = image.width * texCoord[0];
    int j = image.height * texCoord[1];
    return image.pixels[j*image.width + i];
}

glm::vec3 bilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features) {
    glm::vec2 texelPos { image.width * texCoord[0], image.height * texCoord[1]};
    glm::vec2 lowerPos { floor(image.width * texCoord[0]), floor(image.height * texCoord[1]) };
    glm::vec2 upperPos { lowerPos.x + 1, lowerPos.y + 1 };
    float u = texelPos.x - lowerPos.x;
    float v = texelPos.y - lowerPos.y;
    glm::vec3 upperLeft = image.pixels[lowerPos.y * image.width + lowerPos.x];
    glm::vec3 lowerRight = image.pixels[upperPos.y * image.width + upperPos.x];
    glm::vec3 upperRight = image.pixels[upperPos.y * image.width + lowerPos.x];
    glm::vec3 lowerLeft = image.pixels[lowerPos.y * image.width + upperPos.x];
    return upperLeft * (1.0f - u) * (1.0f - v) + lowerRight * u * v + upperRight * (1.0f - u) * v + lowerLeft * u * (1.0f - v);
}