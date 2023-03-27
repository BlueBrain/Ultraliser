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

#include "renderer/Window.h"
#include "renderer/models/NodeTreeModel.h"
#include "renderer/models/VolumeTreeModel.h"

#include <queue>

#include <Ultraliser.h>

class NodeTreeModelFactory
{
public:
    static std::unique_ptr<svorender::Model> create(Ultraliser::SparseOctree &octree)
    {
        auto model = std::make_unique<svorender::NodeTreeModel>();

        auto root = QueueVoxelNode{&octree.getRoot(), octree.getBounds()};
        std::queue<QueueVoxelNode> queue;
        queue.push(root);

        while (!queue.empty())
        {
            auto voxel = queue.front();
            queue.pop();

            if (voxel.node->isLeaf())
            {
                continue;
            }

            auto &min = voxel.bounds.getMin();
            auto &max = voxel.bounds.getMax();
            auto scale = _componentWiseMax(max - min);
            model->addNode(glm::vec3(min.x(), min.y(), min.z()), scale);

            _addNodeChildren(voxel, queue);
        }

        return model;
    }

private:
    struct QueueVoxelNode
    {
        Ultraliser::SparseOctreeNode *node;
        Ultraliser::Bounds bounds;
    };

    static float _componentWiseMax(const Ultraliser::Vector3f &v)
    {
        if (v.x() > v.y() && v.x() > v.z())
        {
            return v.x();
        }
        if (v.y() > v.z())
        {
            return v.y();
        }
        return v.z();
    }

    static void _addNodeChildren(QueueVoxelNode &voxel, std::queue<QueueVoxelNode> &queue)
    {
        auto node = voxel.node;
        for (size_t i = 0; i < node->getNumChildren(); ++i)
        {
            auto &child = node->getChild(i);
            QueueVoxelNode newNode;
            newNode.node = &child;
            newNode.bounds = Ultraliser::SparseOctreeNodeBounds::fromParentBounds(voxel.bounds, child.getSlotMask());
            queue.push(newNode);
        }
    }
};

class VolumeTreeModelFactory
{
public:
    static std::unique_ptr<svorender::Model> create(const Ultraliser::SparseOctree &octree)
    {
        auto &bounds = octree.getBounds();
        auto min = glm::vec3(bounds.getMin().x(), bounds.getMin().y(), bounds.getMin().z());
        auto max = glm::vec3(bounds.getMax().x(), bounds.getMax().y(), bounds.getMax().z());

        auto texture = _generateTexture(octree);

        return std::make_unique<svorender::VolumeTreeModel>(min, max, std::move(texture));
    }

private:
    static svorender::Texture3D _generateTexture(const Ultraliser::SparseOctree &octree)
    {
        auto resolution = static_cast<uint32_t>(1 << octree.getMaxDepth());
        auto values = _sample(octree);
        return svorender::Texture3D(resolution, resolution, resolution, values);
    }

    static std::vector<uint8_t> _sample(const Ultraliser::SparseOctree &tree)
    {
        auto gridResolution = static_cast<uint32_t>(1 << tree.getMaxDepth());
        auto &bounds = tree.getBounds();
        auto srcSize = bounds.getDimensions();
        auto size = glm::vec3(srcSize.x(), srcSize.y(), srcSize.z());
        auto &treeMinBound = bounds.getMin();

        // Cell size
        auto dx = size.x / static_cast<float>(gridResolution);
        auto dy = size.y / static_cast<float>(gridResolution);
        auto dz = size.z / static_cast<float>(gridResolution);

        // Cell center
        auto centerx = dx * 0.5f;
        auto centery = dy * 0.5f;
        auto centerz = dz * 0.5f;

        auto frameLen = gridResolution * gridResolution;

        auto result = std::vector<uint8_t>(gridResolution * gridResolution * gridResolution, 0u);
        auto sampler = Ultraliser::PointSampler(tree);

        for (size_t i = 0; i < gridResolution; i++)
        {
            for (size_t j = 0; j < gridResolution; j++)
            {
                for (size_t k = 0; k < gridResolution; k++)
                {
                    auto localPoint = Ultraliser::Vector3f(dx * k + centerx, dy * j + centery, dz * i + centerz);
                    auto samplePoint = treeMinBound + localPoint;
                    auto filled = sampler.sample(samplePoint);
                    result[frameLen * i + gridResolution * j + k] = filled ? 1 : 0;
                }
            }
        }

        return result;
    }
};

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    auto minS = Ultraliser::Vector3f(-50.0f);
    auto maxS = Ultraliser::Vector3f(50.0f);
    auto bounds = Ultraliser::Bounds(minS, maxS);
    auto depth = static_cast<uint8_t>(8);
    auto octree = Ultraliser::SparseOctree(bounds, depth);
    auto volSizeRaw = bounds.getDimensions();
    auto volumeSize = glm::vec3(volSizeRaw.x(), volSizeRaw.y(), volSizeRaw.z());

    auto pointsToVoxelize = std::vector<Ultraliser::Vector3f>{
        Ultraliser::Vector3f(40.f),
        Ultraliser::Vector3f(-40.f),
        Ultraliser::Vector3f(5.f),
        Ultraliser::Vector3f(100.f),
        Ultraliser::Vector3f(50.0f),
        Ultraliser::Vector3f(-50.0f),
        Ultraliser::Vector3f(0.f)};

    auto pointVoxelizer = Ultraliser::PointVoxelizer(octree);

    for (auto &p : pointsToVoxelize)
    {
        pointVoxelizer.voxelize(p);
    }

    auto a = Ultraliser::Vector3f(0.f, 1.f, 10.f);
    auto b = Ultraliser::Vector3f(5.f, 1.f, 0.f);
    auto c = Ultraliser::Vector3f(-5.f, 1.f, 0.f);

    auto triangleVoxelizer = Ultraliser::TriangleVoxelizer(octree);
    triangleVoxelizer.voxelize(a, b, c);

    // ==================================================================================

    svorender::Window newWindow(1024, 1024);

    newWindow.addModel(NodeTreeModelFactory::create(octree));
    newWindow.addModel(VolumeTreeModelFactory::create(octree));

    newWindow.getCamera().setPosition(glm::vec3(0.f, 0.f, -glm::length(volumeSize)));

    newWindow.renderLoop();

    return 0;
}