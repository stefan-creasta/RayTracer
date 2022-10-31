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

struct Node {
    // Indices in either in meshTrianglePairs or Nodes, the former if the Node is a leaf and the latter if it is not. You can assume children is of size 2 if the node is an internal node.
    std::vector<size_t> children;
    AxisAlignedBox axisAlignedBox;
    bool isLeaf;
    int depth;
    int numLevels;
    int numLeaves;
    float t;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // update parameters when ray hits a triangle
    void triangleIntersectUpdate(const glm::uvec3& tri,  HitInfo& hitInfo, const Ray& ray, const Mesh& mesh, const Features& features) const;

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    // root will be -1 for an empty tree (for example for the "Spheres" scene). Otherwise it will be the index where the root node resides in nodes.
    size_t root;
    std::vector<Node> nodes;
    std::vector<MeshTrianglePair> meshTrianglePairs;
};