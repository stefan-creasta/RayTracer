#pragma once
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <framework/mesh.h>
#include <framework/ray.h>
#include <optional>
#include <variant>
#include <vector>
#include "common.h"
#include <environment_mapping.h>

enum SceneType {
    SingleTriangle,
    Cube,
    CubeTextured,
    CornellBox,
    CornellBoxParallelogramLight,
    Monkey,
    Teapot,
    Dragon,
    Spheres,
    Custom,
    TransparencyDebug,
    TextureDebug
};

struct Scene {
    SceneType type;
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight>> lights;
    std::vector<const EnvironmentMap*> environmentMap;
};

// Load a prebuilt scene.
Scene loadScenePrebuilt(SceneType type, const std::filesystem::path& dataDir);

// Load a scene from a file.
Scene loadSceneFromFile(const std::filesystem::path& path, const std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight>>& lights);
