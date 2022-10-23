#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;
struct MeshTrianglePair {
    Mesh* mesh;
    int triangle;
    glm::vec3 centroid;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Constructor. Used to create a leaf node or an empty BoundingVolumeHierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const std::optional<MeshTrianglePair>& meshTrianglePair, const AxisAlignedBox axisAlignedBox);

    // Constructor. Used to create the inner nodes of the BoundingVolumeHierarchy tree.
    BoundingVolumeHierarchy(Scene* pScene, BoundingVolumeHierarchy& hLeft, BoundingVolumeHierarchy& hRight, const AxisAlignedBox axisAlignedBox);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Return true iff the node is a leaf node.
    [[nodiscard]] bool isLeaf() const;

    // Return if the tree doesn't contain any bounding volumes to traverse.
    [[nodiscard]] bool isEmpty() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);
    void debugDrawLevelHelper(int totalDepth, int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    static BoundingVolumeHierarchy bvhSplitHelper(Scene* pScene, std::vector<BoundingVolumeHierarchy>& allNodes, const std::span<size_t>& indices, const std::span<MeshTrianglePair>& meshTrianglePairs, int direction);


private:
    int m_numLevels;
    int m_numLeaves;
    bool m_isLeaf;
    bool m_isEmpty;
    Scene* m_pScene;
    std::optional<MeshTrianglePair> meshTrianglePair;
    AxisAlignedBox axisAlignedBox;
    std::vector<BoundingVolumeHierarchy> children;
    std::vector<BoundingVolumeHierarchy> allNodes;
};