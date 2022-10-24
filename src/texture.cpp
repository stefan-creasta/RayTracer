#include "texture.h"
#include <framework/image.h>
#include <iostream>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    int i = image.width * texCoord[0];
    int j = image.height * texCoord[1];
    return image.pixels[j*image.width + i];
}