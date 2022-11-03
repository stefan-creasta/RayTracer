#include "sampling.h"
#include <iostream>
#include <random>

// Does multi-jittered sampling of the given axis aligned rectangle.
std::vector<glm::vec2> sample2D(const AxisAlignedRectangle& rect, int h, int k)
{
    std::random_device rd;
    std::mt19937 urng(rd());
    std::uniform_real_distribution<float> distribution(0.f, 1.f);

    std::vector<int> shufflingIndices;
    shufflingIndices.reserve(h);
    for (int i = 0; i < h; i++) { 
        shufflingIndices.push_back(i);
    }
    std::vector<std::vector<glm::vec2>> rows;
    rows.reserve(k);

    for (int j = 0; j < k; j++) {
        std::vector<glm::vec2> row;
        row.reserve(h);
        for (int i = 0; i < h; i++) { 
            row.push_back({j, i});
        }
        rows.push_back(row);
    }

    for (int j = 0; j < k; j++) { 
        std::shuffle(rows[j].begin(), rows[j].end(), urng);
    }

    for (int i = 0; i < h; i++) { 
        std::shuffle(shufflingIndices.begin(), shufflingIndices.end(), urng);

        int j = 0;
        glm::vec2 temp = rows[j][i];
        for (int x = 0; x < k; x++) { 
            rows[j][i] = rows[shufflingIndices[x]][i];
            j = shufflingIndices[x];
        }
        rows[j][i] = temp;
    }

    std::vector<glm::vec2> samples;
    const glm::vec2 scale = glm::vec2 {1. / (k + 1), 1. / (h + 1)};
    const glm::vec2 rectSize = (rect.upper - rect.lower);
    samples.reserve(h * k);
    for (int j = 0; j < k; j++) {
        for (int i = 0; i < h; i++) {
            glm::vec2& center = rows[j][i];
            glm::vec2 jitter = { distribution(urng), distribution(urng) };
            glm::vec2 subdivCoords = center + jitter;
            glm::vec2 unscaledCoords = glm::vec2 { i, j } + scale * subdivCoords;
            samples.push_back(rect.lower + rectSize * (glm::vec2 { i / (h + 1.), j / (k + 1.) } + scale * subdivCoords) );
        }
    }

    return samples;
}
