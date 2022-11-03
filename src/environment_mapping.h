#pragma once
#include <framework/image.h>
#include "common.h"
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <optional>

enum EnvironmentMappingType {
    CYLINDRICAL,
    SPHERICAL,
    SPHEROCYLINDRICAL,
    CUBE
};

class EnvironmentMap {
public:
    static EnvironmentMap loadEnvironmentMap(const std::filesystem::path& file, EnvironmentMappingType mappingType, float verticalFieldOfView, glm::vec3 backgroundColor = glm::vec3 {});
    EnvironmentMap(glm::vec3 color = glm::vec3 {0.f});
    EnvironmentMap(Image& image, EnvironmentMappingType mappingType, float verticalFieldOfView = 180.f, glm::vec3 backgroundColor = glm::vec3 {});

    void buildRadianceHierarchy(const AxisAlignedRectangle& aar, const float totalRadiance);
    std::vector<Ray> getSamplingRay(const glm::vec3& position, const glm::vec3& normal, int n) const;
    Ray getRayForCoordinate(const glm::vec2& coordinates) const;
    glm::vec3 getColor(Ray ray, const Features& features) const;

private:
    const std::optional<Image> texture = std::nullopt;
    const glm::vec3 backgroundColor;
    const float verticalFOVFactor = 1.f;
    const EnvironmentMappingType mappingType = SPHERICAL;
    std::vector<AxisAlignedRectangle> radianceBins;
};