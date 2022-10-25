#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <queue>
#include <iostream>

// Calculate the centroid of a mesh triangle referenced using a MeshTrianglePair.
void calculateCentroid(MeshTrianglePair& meshTrianglePair)
{
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

// Calculate the axis-aligned bounding box of a triangle defined by its vertices.
AxisAlignedBox getTriangleAAB(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2)
{
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

// Calculate the axis-aligned bounding box of a mesh triangle referenced using a MeshTrianglePair.
AxisAlignedBox getTriangleAAB(const MeshTrianglePair& meshTrianglePair)
{
    std::vector<glm::uvec3>& triangles = meshTrianglePair.mesh->triangles;
    std::vector<Vertex>& vertices = meshTrianglePair.mesh->vertices;
    int x = triangles[meshTrianglePair.triangle].x;
    int y = triangles[meshTrianglePair.triangle].y;
    int z = triangles[meshTrianglePair.triangle].z;
    glm::vec3& p0 = vertices[x].position;
    glm::vec3& p1 = vertices[y].position;
    glm::vec3& p2 = vertices[z].position;
    return getTriangleAAB(p0, p1, p2);
}

// Merge two AxisAlignedBoxes into one that contains both with minimal size.
AxisAlignedBox mergeAABs(const AxisAlignedBox& A, const AxisAlignedBox& B)
{
    return AxisAlignedBox {
        glm::min(glm::min(A.lower, B.lower), glm::min(A.upper, B.upper)),
        glm::max(glm::max(A.lower, B.lower), glm::max(A.upper, B.upper)),
    };
}

// Sort the indices according to the triangle centroids in the chosen direction and return the median index.
size_t getMedian(const std::span<size_t>& indices, const std::span<MeshTrianglePair>& meshTrianglePairs, int direction)
{
    std::sort(indices.begin(), indices.end(), [meshTrianglePairs, direction](const size_t a, const size_t b) {
        return (
            (direction == 0) && (meshTrianglePairs[a].centroid.x > meshTrianglePairs[b].centroid.x) || (direction == 1) && (meshTrianglePairs[a].centroid.y > meshTrianglePairs[b].centroid.y) || (direction == 2) && (meshTrianglePairs[a].centroid.z > meshTrianglePairs[b].centroid.z));
    });
    return indices.size() / 2;
}

// Create a Node for the given indices of scene triangles and return its index in nodes.
size_t bvhSplitHelper(Scene* pScene, std::vector<Node>& nodes, const std::span<size_t>& indices, const std::span<MeshTrianglePair>& meshTrianglePairs, int direction, int currentDepth, int maxDepth = -1)
{
    Node result;

    if (maxDepth > -1 && maxDepth <= currentDepth || indices.size() < 2) {
        AxisAlignedBox aab = getTriangleAAB(meshTrianglePairs[indices[0]]);

        for (size_t i : indices) {
            aab = mergeAABs(aab, getTriangleAAB(meshTrianglePairs[i]));
        }

        nodes.push_back(Node { std::vector<size_t>(indices.begin(), indices.end()), aab, true, currentDepth, 1, 1 });
        return nodes.size() - 1;
    }

    size_t median = getMedian(indices, meshTrianglePairs, direction);
    const std::span<size_t> leftIndices = std::span<size_t>(indices.begin(), indices.begin() + median);
    const std::span<size_t> rightIndices = std::span<size_t>(indices.begin() + median, indices.end());

    size_t left = bvhSplitHelper(pScene, nodes, leftIndices, meshTrianglePairs, (direction + 1) % 3, currentDepth + 1, maxDepth);
    size_t right = bvhSplitHelper(pScene, nodes, rightIndices, meshTrianglePairs, (direction + 1) % 3, currentDepth + 1, maxDepth);

    nodes.push_back(Node { 
        { left, right }, 
        mergeAABs(nodes[left].axisAlignedBox, nodes[right].axisAlignedBox), 
        false, 
        currentDepth,
        glm::max(nodes[left].numLevels, nodes[right].numLevels) + 1, 
        nodes[left].numLeaves + nodes[right].numLeaves 
    });
    return nodes.size() - 1;
}

// Constructor. Receives the scene and builds the bounding volume hierarchy.
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    size_t n = 0;
    for (Mesh& mesh : pScene->meshes) {
        n += mesh.triangles.size();
    }

    nodes = std::vector<Node>();

    if (n == 0) {
        root = -1;
        m_numLevels = 0;
        m_numLeaves = 0;
        return;
    }

    std::vector<size_t> indices = std::vector<size_t>();
    meshTrianglePairs = std::vector<MeshTrianglePair>();
    indices.reserve(n);
    meshTrianglePairs.reserve(n);
    nodes.reserve(2 * n - 1);

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

    root = bvhSplitHelper(pScene, nodes, indices, meshTrianglePairs, 0, 0);
    m_numLevels = nodes[root].numLevels;
    m_numLeaves = nodes[root].numLeaves;
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

void debugDrawLevelHelper(std::vector<Node>& nodes, size_t node, int level)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    if (nodes[node].depth == level) {
        drawAABB(nodes[node].axisAlignedBox, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (!nodes[node].isLeaf) {
            // Adjust total depth by how much each subtree differs in height from the one with greater height.
            debugDrawLevelHelper(nodes, nodes[node].children[0], level);
            debugDrawLevelHelper(nodes, nodes[node].children[1], level);
        }
    }
    // drawAABB(aabb, DrawMode::Wireframe);
}
// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    if (root == -1)
        return;

    debugDrawLevelHelper(nodes, root, level);
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    if (root == -1)
        return;

    Node current = nodes[root];
    int currentLeafIdx = leafIdx;

    while (!current.isLeaf) { 
        if (nodes[current.children[0]].numLeaves > currentLeafIdx) {
            current = nodes[current.children[0]];
        } else {
            currentLeafIdx -= nodes[current.children[0]].numLeaves;
            current = nodes[current.children[1]];
        }
    }

    drawAABB(current.axisAlignedBox, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    // drawAABB(aabb, DrawMode::Wireframe);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
/** struct myComp {
    constexpr bool operator()(
        Node const& a,
        Node const& b)
        const noexcept
    {
        return a.t < b.t;
    }
};**/
float INF = std::numeric_limits<float>::infinity();
class Compare {
public:
    bool operator()(Node a, Node b) {
        return a.t > b.t;
    }
};
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    bool hit = false;
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres)
        hit |= intersectRayWithShape(sphere, ray, hitInfo);
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
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
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        std::priority_queue<Node, std::vector<Node>, Compare> pq;
        if (root == -1) {
            return hit;
        }
        pq.push(nodes[root]);
        bool hitTri = false;
        float sphereT = ray.t;
        float minT = INF;
        int minTri = -1;
        while (!pq.empty()) {
            Node front = pq.top();
            pq.pop();
            if (minT < front.t) {
                drawAABB(front.axisAlignedBox, DrawMode::Wireframe, glm::vec3 {1.0f, 0.0f, 0.0f});
                continue;
            }
            drawAABB(front.axisAlignedBox, DrawMode::Wireframe);
            if (front.isLeaf == true) {
                ray.t = INF;
                for (size_t currentChild : front.children) {
                    MeshTrianglePair pair = meshTrianglePairs[currentChild];
                    Mesh mesh = *pair.mesh;
                    glm::uvec3 tri = mesh.triangles[pair.triangle];
                    const auto v0 = mesh.vertices[tri[0]];
                    const auto v1 = mesh.vertices[tri[1]];
                    const auto v2 = mesh.vertices[tri[2]];
                    if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                        hitInfo.material = mesh.material;
                        hitTri = true;
                        if (ray.t < minT)
                            minTri = currentChild;
                    }
                }
                if (hitTri == true) {
                    if (ray.t < minT) {
                        minT = ray.t;
                    }
                }
            } else {
                Node left = nodes[front.children[0]];
                Node right = nodes[front.children[1]];
                ray.t = INF;
                if (intersectRayWithShape(left.axisAlignedBox, ray) == true) {
                    left.t = ray.t;
                    pq.push(left);
                }
                ray.t = INF;
                if (intersectRayWithShape(right.axisAlignedBox, ray) == true) {
                    right.t = ray.t;
                    pq.push(right);
                }
            }
        }
        if (ray.t > minT) {
            ray.t = minT;
        }
        if (sphereT > minT) {
            MeshTrianglePair pair = meshTrianglePairs[minTri];
            Mesh mesh = *pair.mesh;
            glm::uvec3 tri = mesh.triangles[pair.triangle];
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            drawTriangle(v0, v1, v2);
        }
        return hit | hitTri;
    }
    return hit;
}