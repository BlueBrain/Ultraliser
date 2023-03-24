/* Copyright (c) 2020, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
 *
 * This file is part of SimCrusher
 * <LINK>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <iostream>

#include "renderer/Common.h"
#include "renderer/Window.h"
#include "renderer/materials/BoundRenderMaterial.h"
#include "renderer/materials/VolumeMaterial.h"

#include <cstdlib>
#include <queue>

#include <Ultraliser.h>

class SparseNodeBounds
{
public:
    static Ultraliser::Bounds fromParentBounds(const Ultraliser::Bounds &parentBounds, uint8_t slotMask)
    {
        auto left = slotMask & Ultraliser::SparseOctreeNodeSlot::backBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::backTopLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontTopLeft;
        auto bottom = slotMask & Ultraliser::SparseOctreeNodeSlot::backBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::backBottomRight
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontBottomRight;
        auto back = slotMask < Ultraliser::SparseOctreeNodeSlot::frontBottomLeft;

        auto xMult = left ? 0.f : 1.f;
        auto yMult = bottom ? 0.f : 1.f;
        auto zMult = back ? 0.f : 1.f;
        auto mult = Ultraliser::Vector3f(xMult, yMult, zMult);

        auto min = parentBounds.getMin();
        auto center = parentBounds.getCenter();
        auto halfSize = center - min;

        auto childMin = min + halfSize * mult;
        auto childMax = center + halfSize * mult;
        return Ultraliser::Bounds(childMin, childMax);
    }
};

struct QueueVoxelNode
{
    Ultraliser::SparseOctreeNode *node;
    Ultraliser::Bounds bounds;
};

std::vector<scr::CubeMesh> getIntermediateVoxelMehes(Ultraliser::SparseOctree &octree)
{
    Ultraliser::SparseOctreeNode *node = &octree.getRoot();
    QueueVoxelNode root;
    root.node = node;
    root.bounds = octree.getBounds();
    std::queue<QueueVoxelNode> q;
    q.push(root);

    std::vector<scr::CubeMesh> result;

    while (!q.empty())
    {
        QueueVoxelNode n = q.front();
        q.pop();

        // Only add if its not a leaf node
        if (n.node->isLeaf())
        {
            continue;
        }

        auto &min = n.bounds.getMin();
        auto &max = n.bounds.getMax();
        result.emplace_back(glm::vec3(min.x(), min.y(), min.z()), glm::vec3(max.x(), max.y(), max.z()));

        for (size_t i = 0; i < n.node->getNumChildren(); ++i)
        {
            auto &childNode = n.node->getChild(i);
            QueueVoxelNode newQueueNode;
            newQueueNode.node = &childNode;
            newQueueNode.bounds = SparseNodeBounds::fromParentBounds(n.bounds, childNode.getSlotMask());
            q.push(newQueueNode);
        }
    }

    return result;
}

// Samples the octree into a regular grid (3D Texture)
std::vector<uint8_t> octreeToVolumeData(const sc::SparseOctree &tree)
{
    const uint32_t gridResolution = 1 << tree.getMaxDepth();
    const sc::Point3DF &srcSize = tree.get3DSize();
    const glm::vec3 size(srcSize.x, srcSize.y, srcSize.z);

    const sc::Point3DF treeMinBound = tree.getBounds().min;

    // Cell size
    const float dx = size.x / static_cast<float>(gridResolution);
    const float dy = size.y / static_cast<float>(gridResolution);
    const float dz = size.z / static_cast<float>(gridResolution);

    // Cell center
    const float centerx = dx / 2.f;
    const float centery = dy / 2.f;
    const float centerz = dz / 2.f;

    const uint32_t frameLen = gridResolution * gridResolution;

    const sc::VoxelPointTest algorithm;

    std::vector<uint8_t> result;
    result.resize(gridResolution * gridResolution * gridResolution, 0u);

    for (size_t i = 0; i < gridResolution; i++)
    {
        for (size_t j = 0; j < gridResolution; j++)
        {
            for (size_t k = 0; k < gridResolution; k++)
            {
                const sc::Point3DF samplePoint =
                    treeMinBound + sc::Point3DF(dx * k + centerx, dy * j + centery, dz * i + centerz);

                const uint8_t sample = sc::sampleShape(tree, algorithm, sc::PointShape(samplePoint));

                result[frameLen * i + gridResolution * j + k] = sample;
            }
        }
    }

    return result;
}

int main(SCR_UNUSED int argc, SCR_UNUSED char **argv)
{
    std::cout << "Creating Octree..." << std::endl;

    const sc::Point3DF minS(-50.0f);
    const sc::Point3DF maxS(50.0f);
    const uint8_t depth = 8;
    sc::SparseOctree octree(minS, maxS, depth);
    const sc::Point3DF volSizeRaw = octree.get3DSize();
    const glm::vec3 volumeSize(volSizeRaw.x, volSizeRaw.y, volSizeRaw.z);

    // Points to voxelize
    std::vector<sc::Point3DF> pointsToSample = {
        sc::Point3DF(40.f),
        sc::Point3DF(-40.f),
        sc::Point3DF(5.f),
        sc::Point3DF(100.f),
        sc::Point3DF(50.0f),
        sc::Point3DF(-50.0f),
        sc::Point3DF(0.f)};

    sc::VoxelPointTest algorithm;
    for (const auto &p : pointsToSample)
        sc::voxelizeShape(octree, algorithm, sc::PointShape(p), rand() % 255);

    const float minBound = 0.0f;
    const float maxBound = 10.0f;
    const float depth = 8;
    const auto minS = sc::Point3DF(minBound);
    const auto maxS = sc::Point3DF(maxBound);
    const float voxelSize = 0.0390625f;
    sc::SparseOctree octree(minS, maxS, depth);
    sc::VoxelPointTest voxelPointTest;
    const uint8_t testValue = 128u;
    const glm::vec3 volumeSize(maxBound - minBound);

    sc::voxelizeShape(octree, voxelPointTest, sc::PointShape(sc::Point3DF(minBound + voxelSize * 0.5f)), testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 1.5f, minBound + voxelSize * 0.5f, minBound + voxelSize * 0.5f)),
        testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 1.5f, minBound + voxelSize * 1.5f, minBound + voxelSize * 0.5f)),
        testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 0.5f, minBound + voxelSize * 1.5f, minBound + voxelSize * 0.5f)),
        testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 0.5f, minBound + voxelSize * 0.5f, minBound + voxelSize * 1.5f)),
        testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 1.5f, minBound + voxelSize * 0.5f, minBound + voxelSize * 1.5f)),
        testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 1.5f, minBound + voxelSize * 1.5f, minBound + voxelSize * 1.5f)),
        testValue);
    sc::voxelizeShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 0.5f, minBound + voxelSize * 1.5f, minBound + voxelSize * 1.5f)),
        0u); // should be zero

    octree.compact();

    auto sample = sc::sampleShape(
        octree,
        voxelPointTest,
        sc::PointShape(
            sc::Point3DF(minBound + voxelSize * 0.5f, minBound + voxelSize * 1.5f, minBound + voxelSize * 1.5f)));
    std::cout << "Sampled: " << static_cast<int>(sample) << std::endl;

    std::cout << "Sparse Octree:" << std::endl;
    octree.printStats();

    sc::GPUSparseOctree gpuOctree(octree);
    std::cout << "GPU Sparse Octree:" << std::endl;
    gpuOctree.printStats();

    // ==================================================================================

    scr::Window newWindow(1024, 1024);

    // BOUNDING BOXES FOR INTERMEDIATE VOXELS
    scr::ShaderConfig wire;
    wire.vertexShader = "./shaders/boundrenderV.glsl";
    wire.geometryShader = "./shaders/boundrenderG.glsl";
    wire.fragmentShader = "./shaders/boundrenderF.glsl";
    scr::Shader wireShader(wire);
    scr::BoundRenderMaterial wireMat(&wireShader);
    std::vector<scr::CubeMesh> parentVoxelMeshes = getIntermediateVoxelMehes(octree);
    for (auto &parentMeshes : parentVoxelMeshes)
    {
        newWindow.addModel(&parentMeshes, &wireMat);
    }

    // VOLUME
    scr::CubeMesh boundBox(glm::vec3(minS.x, minS.y, minS.z), glm::vec3(maxS.x, maxS.y, maxS.z));
    scr::ShaderConfig sc;
    sc.vertexShader = "./shaders/volrenderV.glsl";
    sc.fragmentShader = "./shaders/volrenderF.glsl";
    scr::Shader shader(sc);
    const uint32_t resolution = 1 << octree.getMaxDepth();
    const std::vector<uint8_t> volData = octreeToVolumeData(octree);
    scr::Texture3D volume(resolution, resolution, resolution, volData);
    scr::VolumeMaterial volumeMat(&shader, &volume);
    volumeMat.setBounds(glm::vec3(minS.x, minS.y, minS.z), glm::vec3(maxS.x, maxS.y, maxS.z));
    newWindow.addModel(&boundBox, &volumeMat);

    // CAMERA
    const float len = glm::length(volumeSize);
    newWindow.getCamera().setPosition(glm::vec3(0.f, 0.f, -len));

    newWindow.addModel(&boundBox, &wireMat);

    newWindow.renderLoop();

    return 0;
}
* /

    int main()
{
    return 0;
}