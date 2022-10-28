#include "texture.h"
#include <framework/image.h>
#include <iostream>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    int i = (image.width * texCoord[0] - 0.5f);
    int j = (image.height * texCoord[1] - 0.5f);

    // clamp if i and j are out of range

    i = std::max(0, std::min(image.width - 1, i));
    j = std::max(0, std::min(image.height - 1, j));
    return image.pixels[j*image.width + i];
}