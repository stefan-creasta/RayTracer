#include "texture.h"
#include "draw.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <iostream>

ImageMipMap getMipMap(const Image& image)
{
    ImageMipMap mipmap;
    mipmap.height.push_back(image.height);
    mipmap.width.push_back(image.width);
    mipmap.pixels.push_back(image.pixels);
    int last = 0;
    while (1) {
        if (mipmap.pixels[last].size() <= 1) {
            break;
        }
        std::vector<glm::vec3> vec;
        //std::cout << mipmap.height[last] << std::endl;
        for (int i = 0; i < mipmap.height[last]; i += 2) {
            for (int j = 0; j < mipmap.width[last]; j += 2) {
                image.pixels[j * image.width + i];
                glm::vec3 avg = mipmap.pixels[last][i * mipmap.width[last] + j] + mipmap.pixels[last][i * mipmap.width[last] + j + 1] + mipmap.pixels[last][(i + 1) * mipmap.width[last] + j] + mipmap.pixels[last][(i + 1) * mipmap.width[last] + j + 1];
                avg *= (1.0f / 4.0f);
                vec.push_back(avg);
            }
        }
        mipmap.height.push_back(mipmap.height[last] / 2);
        mipmap.width.push_back(mipmap.width[last] / 2);
        mipmap.pixels.push_back(vec);
        last++;
    }
    return mipmap;
}

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 position = ray.origin + ray.t * ray.direction;
    glm::vec3 lightDir = glm::normalize(lightPosition - position);
    float dot = glm::dot(glm::normalize(hitInfo.normal), lightDir);
    float eps = 1e-6;
    if (dot < eps) {
        return {0, 0, 0};
    }
    glm::vec3 view = ray.origin - position;
    glm::vec3 h = glm::normalize(lightDir + view);
    if (glm::dot(h, glm::normalize(hitInfo.normal)) < eps) {
        return {0, 0, 0};
    }
    // EXPERIMENT
    if (hitInfo.material.kdTexture && features.enableTextureMapping) {
        glm::vec3 texel = acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
        //ImageMipMap mipmap = getMipMap(*hitInfo.material.kdTexture.get());
        if (features.extra.enableBilinearTextureFiltering) {

            texel = bilinearInterpolation(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
            glm::vec2 debugUV = getUVForBilinear(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
            glm::vec3 w = glm::normalize(hitInfo.normal);
            glm::vec3 t = glm::normalize(w - glm::vec3 { 0.1f, 0.0f, 0.0f });
            glm::vec3 u = glm::normalize(glm::cross(t, w));
            glm::vec3 v = glm::normalize(glm::cross(w, u));
            Ray rayU = Ray { position, u, debugUV.x / 5.0f };
            Ray rayV = Ray { position, v, debugUV.y / 5.0f };
            Ray rayOU = Ray { position, -u, (1 - debugUV.x) / 5.0f };
            Ray rayOV = Ray { position, -v, (1 - debugUV.y) / 5.0f };
            drawRay(rayU, glm::vec3 { 1.0f, 0.5f, 0.0f });
            drawRay(rayV, glm::vec3 { 1.0f, 0.0f, 0.5f });
            drawRay(rayOU, glm::vec3 { 1.0f, 0.5f, 0.0f });
            drawRay(rayOV, glm::vec3 { 1.0f, 0.0f, 0.5f });
            //std::cout << rayU.origin.x << std::endl;
        }
        return lightColor * hitInfo.material.ks * pow(glm::dot(glm::normalize(hitInfo.normal), h), hitInfo.material.shininess) + lightColor * texel * dot; 
    }
    return lightColor * hitInfo.material.ks * pow(glm::dot(glm::normalize(hitInfo.normal), h), hitInfo.material.shininess) + lightColor * hitInfo.material.kd * dot;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    glm::vec3 point = ray.origin + ray.direction * ray.t;
    //glm::vec3 r = glm::reflect(glm::normalize(ray.direction), glm::normalize(hitInfo.normal));
    glm::vec3 r = glm::normalize(ray.direction) - 2 * glm::dot(glm::normalize(hitInfo.normal), glm::normalize(ray.direction)) * glm::normalize(hitInfo.normal);
    Ray reflectionRay {point + float(1e-5)*r, r};
    //drawRay(Ray{point, hitInfo.normal, 1}, glm::vec3{0,1,0.5});
    return reflectionRay;
}