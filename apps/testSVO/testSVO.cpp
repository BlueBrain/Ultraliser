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
            auto &voxel = queue.front();

            auto &min = voxel.bounds.getMin();
            auto &max = voxel.bounds.getMax();
            auto scale = _componentWiseMax(max - min);
            model->addNode(glm::vec3(min.x(), min.y(), min.z()), scale);
            _addNodeChildren(voxel, queue);

            queue.pop();
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
        auto sampler = Ultraliser::OctreePointSampler(tree);

#pragma omp parallel for
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

        result[0] = 1;

        return result;
    }
};

class MorphologyOctreeFactory
{
public:
    static Ultraliser::SparseOctree create(uint8_t maxDepth)
    {
        auto morphology = _readMorphology();
        auto bounds = _getMorhpologyBounds(*morphology);
        auto resolution = 1u << maxDepth;
        auto octree = Ultraliser::SparseOctree(bounds, maxDepth);

        Ultraliser::OctreeNeuronMorphologyVoxelizer::voxelize(octree, *morphology);

        return octree;
    }

private:
    static std::unique_ptr<Ultraliser::NeuronMorphology> _readMorphology()
    {
        auto morphologyPath = "/home/nadir/Desktop/dend-oh140807_A0_idG_axon-mpg141017_a1-2_idC.swc";
        auto reader = Ultraliser::NeuronSWCReader(morphologyPath);
        return std::unique_ptr<Ultraliser::NeuronMorphology>(reader.getMorphology());
    }

    static Ultraliser::Bounds _getMorhpologyBounds(const Ultraliser::NeuronMorphology &morphology)
    {
        auto min = Ultraliser::Vector3f();
        auto max = Ultraliser::Vector3f();
        auto dim = Ultraliser::Vector3f();
        auto center = Ultraliser::Vector3f();
        morphology.getBoundingBox(min, max, dim, center);

        size_t greatest = 0;
        for (size_t i = 1; i < 3; ++i)
        {
            if (dim[i] > dim[greatest])
            {
                greatest = i;
            }
        }

        for (size_t i = 0; i < 3; ++i)
        {
            auto halfDiff = (dim[greatest] - dim[i]) * 0.5f;
            min[i] -= halfDiff;
            max[i] += halfDiff;
        }

        return Ultraliser::Bounds(min, max);
    }
};

class MeshOctreeFactory
{
public:
    static Ultraliser::SparseOctree create(uint8_t depth)
    {
        auto path = "/home/nadir/Desktop/morphology_1.obj";
        auto mesh = Ultraliser::Mesh(path);

        auto min = Ultraliser::Vector3f();
        auto max = Ultraliser::Vector3f();
        mesh.computeBoundingBox(min, max);
        auto bounds = Ultraliser::Bounds(min, max);

        auto octree = Ultraliser::SparseOctree(bounds, depth);

        Ultraliser::OctreeMeshVoxelizer::voxelize(octree, mesh);

        return octree;
    }
};

class TestCompactOctreeFactory
{
public:
    static Ultraliser::SparseOctree create()
    {
        auto min = Ultraliser::Vector3f(-40.f);
        auto max = Ultraliser::Vector3f(40.f);
        auto bounds = Ultraliser::Bounds(min, max);
        auto depth = static_cast<uint8_t>(3);
        auto octree = Ultraliser::SparseOctree(bounds, depth);

        auto voxelizer = Ultraliser::OctreePointVoxelizer(octree);

        auto resolution = 1 << depth;
        auto start = Ultraliser::Vec3ui_32(4);
        auto end = Ultraliser::Vec3ui_32(resolution);

        for (uint32_t x = start.x(); x < end.x(); ++x)
        {
            for (uint32_t y = start.y(); y < end.y(); ++y)
            {
                for (uint32_t z = start.z(); z < end.z(); ++z)
                {
                    auto point = Ultraliser::Vec3ui_32(x, y, z);
                    voxelizer.voxelize(point);
                }
            }
        }

        return octree;
    }
};

class TestFillableOctreeFactory
{
public:
    static Ultraliser::SparseOctree create()
    {
        auto min = Ultraliser::Vector3f(-40.f);
        auto max = Ultraliser::Vector3f(40.f);
        auto bounds = Ultraliser::Bounds(min, max);
        auto depth = static_cast<uint8_t>(3);
        auto octree = Ultraliser::SparseOctree(bounds, depth);

        auto voxelizer = Ultraliser::OctreePointVoxelizer(octree);

        auto resolution = 1 << depth;
        auto start = Ultraliser::Vec3ui_32(4);
        auto end = Ultraliser::Vec3ui_32(resolution);

        for (uint32_t x = start.x(); x < end.x(); ++x)
        {
            for (uint32_t y = start.y(); y < end.y(); ++y)
            {
                for (uint32_t z = start.z(); z < end.z(); ++z)
                {
                    auto borderX = (x == start.x() || x == end.x() - 1);
                    auto borderY = (y == start.y() || y == end.y() - 1);
                    auto borderZ = (z == start.x() || z == end.z() - 1);
                    if (!borderX && !borderY && !borderZ)
                    {
                        continue;
                    }

                    auto point = Ultraliser::Vec3ui_32(x, y, z);
                    voxelizer.voxelize(point);
                }
            }
        }

        return octree;
    }
};

class TestOctreeFactory
{
public:
    static Ultraliser::SparseOctree create()
    {
        auto min = Ultraliser::Vector3f(-40.f);
        auto max = Ultraliser::Vector3f(40.f);
        auto bounds = Ultraliser::Bounds(min, max);
        auto depth = static_cast<uint8_t>(3);
        auto octree = Ultraliser::SparseOctree(bounds, depth);

        auto voxelizer = Ultraliser::OctreePointVoxelizer(octree);

        auto resolution = 1 << depth;

        auto testPoint = Ultraliser::Vector3f(30.f);
        voxelizer.voxelize(testPoint);
        voxelizer.voxelize(testPoint);

        return octree;
    }
};

class OctreeOptimizer
{
public:
    static void optimize(Ultraliser::SparseOctree &octree)
    {
        std::cout << "Before compacting" << std::endl;
        Ultraliser::OctreeStatsPrinter::print(octree);

        std::cout << "After compacting" << std::endl;
        Ultraliser::OctreeCompacter::compact(octree);
        Ultraliser::OctreeStatsPrinter::print(octree);

        std::cout << "After filling" << std::endl;
        Ultraliser::OctreeFloodFiller::fill(octree);
        Ultraliser::OctreeStatsPrinter::print(octree);

        std::cout << "After compacting" << std::endl;
        Ultraliser::OctreeCompacter::compact(octree);
        Ultraliser::OctreeStatsPrinter::print(octree);
    }
};

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    auto octree = MorphologyOctreeFactory::create(10);
    OctreeOptimizer::optimize(octree);

    svorender::Window newWindow(1024, 1024);

    newWindow.addModel(NodeTreeModelFactory::create(octree));
    newWindow.addModel(VolumeTreeModelFactory::create(octree));

    auto volSizeRaw = octree.getBounds().getDimensions();
    auto volumeSize = glm::vec3(volSizeRaw.x(), volSizeRaw.y(), volSizeRaw.z());
    newWindow.getCamera().setPosition(glm::vec3(0.f, 0.f, -glm::length(volumeSize)));

    newWindow.renderLoop();

    return 0;
}