#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>

BoundingVolumeHierarchy* lastBVH;

void calculateCentroid(MeshTrianglePair& meshTrianglePair) { 
    std::vector<glm::uvec3>& triangles = meshTrianglePair.mesh->triangles;
    std::vector<Vertex>& vertices = meshTrianglePair.mesh->vertices;
    int x = triangles[meshTrianglePair.triangle].x;
    int y = triangles[meshTrianglePair.triangle].y;
    int z = triangles[meshTrianglePair.triangle].z;
    glm::vec3& p0 = vertices[x].position;
    glm::vec3& p1 = vertices[y].position;
    glm::vec3& p2 = vertices[z].position;
    meshTrianglePair.centroid = (1.f / 3.f) * (p0 + p1 + p2);
}

AxisAlignedBox getTriangleAAB(const glm::vec3 p0, const glm::vec3 p1, const glm::vec3 p2) { 
    return AxisAlignedBox {
        glm::vec3 {
            std::min({ p0.x, p1.x, p2.x }),
            std::min({ p0.y, p1.y, p2.y }),
            std::min({ p0.z, p1.z, p2.z }),
        },
        glm::vec3 {
            std::max({ p0.x, p1.x, p2.x }),
            std::max({ p0.y, p1.y, p2.y }),
            std::max({ p0.z, p1.z, p2.z }),
        }
    };
}

AxisAlignedBox getTriangleAAB(const MeshTrianglePair& meshTrianglePair) { 
    std::vector<glm::uvec3>& triangles = meshTrianglePair.mesh->triangles;
    std::vector<Vertex>& vertices = meshTrianglePair.mesh->vertices;
    int x = triangles[meshTrianglePair.triangle].x;
    int y = triangles[meshTrianglePair.triangle].y;
    int z = triangles[meshTrianglePair.triangle].z;
    glm::vec3& p0 = vertices[x].position;
    glm::vec3& p1 = vertices[y].position;
    glm::vec3& p2 = vertices[z].position;
    return AxisAlignedBox
    {
        glm::vec3 { 
            std::min({ p0.x, p1.x, p2.x }),
            std::min({ p0.y, p1.y, p2.y }),
            std::min({ p0.z, p1.z, p2.z }),
        },
        glm::vec3 {
            std::max({ p0.x, p1.x, p2.x }),
            std::max({ p0.y, p1.y, p2.y }),
            std::max({ p0.z, p1.z, p2.z }),
        }
    };
}

AxisAlignedBox mergeAABs(const AxisAlignedBox& A, const AxisAlignedBox& B) { 
    return AxisAlignedBox {
        glm::min(glm::min(A.lower, B.lower), glm::min(A.upper, B.upper)),
        glm::max(glm::max(A.lower, B.lower), glm::max(A.upper, B.upper)),
    };
}

size_t getMedian(const std::span<size_t>& indices, const std::span<MeshTrianglePair>& meshTrianglePairs, int direction) { 
    std::sort(indices.begin(), indices.end(), [meshTrianglePairs, direction](const size_t a, const size_t b) {
        return (
            (direction == 0) && (meshTrianglePairs[a].centroid.x > meshTrianglePairs[b].centroid.x) ||
            (direction == 1) && (meshTrianglePairs[a].centroid.y > meshTrianglePairs[b].centroid.y) ||
            (direction == 2) && (meshTrianglePairs[a].centroid.z > meshTrianglePairs[b].centroid.z)
        );
    });
    return indices.size() / 2;
}

BoundingVolumeHierarchy* BoundingVolumeHierarchy::bvhSplitHelper(Scene* pScene, const std::span<size_t>& indices, const std::span<MeshTrianglePair>& meshTrianglePairs, int direction)
{
    if (indices.size() == 0) {
        BoundingVolumeHierarchy* current = new BoundingVolumeHierarchy(pScene, std::optional<MeshTrianglePair>(), AxisAlignedBox {});
        current->m_isEmpty = true;
        return current;
    }
    if (indices.size() == 1) {
        MeshTrianglePair pair = meshTrianglePairs[indices[0]];
        const auto aab = getTriangleAAB(pair);
        BoundingVolumeHierarchy* current = new BoundingVolumeHierarchy(pScene, pair, aab);
        current->m_isEmpty = false;
        return current;
    }
    size_t median = getMedian(indices, meshTrianglePairs, direction);
    const std::span<size_t> leftIndices = std::span<size_t>(indices.begin(), indices.begin() + median);
    const std::span<size_t> rightIndices = std::span<size_t>(indices.begin() + median, indices.end());

    BoundingVolumeHierarchy* hLeft = bvhSplitHelper(pScene, leftIndices, meshTrianglePairs, (direction + 1) % 3);
    BoundingVolumeHierarchy* hRight = bvhSplitHelper(pScene, rightIndices, meshTrianglePairs, (direction + 1) % 3);

    AxisAlignedBox currentBox = mergeAABs(hLeft->axisAlignedBox, hRight->axisAlignedBox);

    BoundingVolumeHierarchy* current = new BoundingVolumeHierarchy(pScene, hLeft, hRight, currentBox);
    current->m_isEmpty = false;
    
    return current;
}

// Constructor. Used to create a leaf node or an empty BoundingVolumeHierarchy.
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const std::optional<MeshTrianglePair>& meshTrianglePair, const AxisAlignedBox axisAlignedBox)
    : m_pScene(pScene)
    , meshTrianglePair(meshTrianglePair)
    , axisAlignedBox(axisAlignedBox)
    , m_numLeaves(1)
    , m_numLevels(1)
    , m_isLeaf(true)
    , m_isEmpty(false)
    , leftChild(nullptr)
    , rightChild(nullptr)
{ 
}

// Constructor. Used to create the inner nodes of the BoundingVolumeHierarchy tree.
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, BoundingVolumeHierarchy* hLeft, BoundingVolumeHierarchy* hRight, const AxisAlignedBox axisAlignedBox)
    : m_pScene(pScene)
    , meshTrianglePair(meshTrianglePair)
    , axisAlignedBox(axisAlignedBox)
    , m_isLeaf(false)
    , m_isEmpty(false)
    , leftChild(hLeft)
    , rightChild(hRight)
{
    this->m_numLeaves = hLeft->m_numLeaves + hRight->m_numLeaves;
    this->m_numLevels = glm::max(hLeft->m_numLevels, hRight->m_numLevels) + 1;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    size_t n = 0;
    for (Mesh& mesh : pScene->meshes) { 
        n += mesh.triangles.size();
    }

    if (n == 0) {
        new (this) BoundingVolumeHierarchy(pScene, std::optional<MeshTrianglePair>(), AxisAlignedBox {});
        this->m_isEmpty = true;
        return;
    }

    std::vector<size_t> indices = std::vector<size_t>();
    std::vector<MeshTrianglePair> meshTrianglePairs = std::vector<MeshTrianglePair>();
    indices.reserve(n);
    meshTrianglePairs.reserve(n);

    // i is incremented in the inner loop!
    size_t i = 0;
    for (Mesh& mesh : pScene->meshes) {
        size_t len = mesh.triangles.size();
        for (int j = 0; j < len; j++) {
            indices.push_back(i);
            MeshTrianglePair pair = MeshTrianglePair { &mesh, j, glm::vec3 {} };
            calculateCentroid(pair);
            meshTrianglePairs.push_back(pair);
            i++;
        }
    }

    if (lastBVH)
        delete lastBVH;

    // change object pointed at by this to our properly created BVH
    *this = *bvhSplitHelper(pScene, indices, meshTrianglePairs, 0);
    lastBVH = this;
}

// Destructor for BoundingVolumeHierarchy
BoundingVolumeHierarchy::~BoundingVolumeHierarchy()
{
    if (this->leftChild)
        delete this->leftChild;
    if (this->rightChild)
        delete this->rightChild;
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return this->m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return this->m_numLeaves;
}

// Return if the node is a leaf.
bool BoundingVolumeHierarchy::isLeaf() const
{
    return this->m_isLeaf;
}

// Return if the tree doesn't contain any bounding volumes to traverse.
bool BoundingVolumeHierarchy::isEmpty() const
{
    return this->m_isEmpty;
}

void BoundingVolumeHierarchy::debugDrawLevelHelper(int totalDepth, int level) {
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    if (totalDepth - this->numLevels() >= level) {
        drawAABB(axisAlignedBox, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (!this->isLeaf()) {
            // Adjust total depth by how much each subtree differs in height from the one with greater height.
            leftChild->debugDrawLevelHelper(totalDepth - (this->numLevels() - 1 - leftChild->numLevels()), level);
            rightChild->debugDrawLevelHelper(totalDepth - (this->numLevels() - 1 - rightChild->numLevels()), level);
        }
    }
    // drawAABB(aabb, DrawMode::Wireframe);
}
// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    debugDrawLevelHelper(this->numLevels(), level);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    if (this->isLeaf()) {
        drawAABB(axisAlignedBox, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (leftChild->numLeaves() > leafIdx)
            leftChild->debugDrawLeaf(leafIdx);
        else
            rightChild->debugDrawLeaf(leafIdx - leftChild->numLeaves());
    }
    //drawAABB(aabb, DrawMode::Wireframe);
    

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}