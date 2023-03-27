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

#include "NodeTreeModel.h"

namespace
{
class ShaderFactory
{
public:
    static svorender::Shader create()
    {
        svorender::ShaderConfig wire;
        wire.vertexShader = "./shaders/boundrenderV.glsl";
        wire.geometryShader = "./shaders/boundrenderG.glsl";
        wire.fragmentShader = "./shaders/boundrenderF.glsl";
        return svorender::Shader(wire);
    }
};

class NodeDataBufferFactory
{
public:
    static void create(GLuint &ssb, const std::vector<glm::vec4> &nodes)
    {
        if (ssb != GL_INVALID_VALUE)
        {
            glDeleteBuffers(1, &ssb);
        }

        auto size = nodes.size() * sizeof(glm::vec4);
        glGenBuffers(1, &ssb);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssb);
        glBufferData(GL_SHADER_STORAGE_BUFFER, size, nodes.data(), GL_STATIC_DRAW);
    }
};

struct DrawCommand
{
    uint32_t count;
    uint32_t instanceCount;
    uint32_t firstIndex = 0;
    int baseVertex = 0;
    uint32_t baseInstance = 0;
};

class DrawCommandBufferFactory
{
public:
    static void create(GLuint &buffer, const DrawCommand &command)
    {
        if (buffer != GL_INVALID_VALUE)
        {
            glDeleteBuffers(1, &buffer);
        }

        glGenBuffers(1, &buffer);
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, buffer);
        glBufferData(GL_DRAW_INDIRECT_BUFFER, sizeof(DrawCommand), &command, GL_STATIC_READ);
    }
};
}

namespace svorender
{
NodeTreeModel::NodeTreeModel()
    : _shader(ShaderFactory::create())
{
}

void NodeTreeModel::addNode(const glm::vec3 &position, float scale)
{
    _dirty = true;
    _nodes.emplace_back(position, scale);
}

void NodeTreeModel::render(const Camera &camera)
{
    if (_dirty)
    {
        _dirty = false;
        NodeDataBufferFactory::create(_nodeDataBuffer, _nodes);

        auto command = DrawCommand{36, static_cast<uint32_t>(_nodes.size())};
        DrawCommandBufferFactory::create(_drawCommandBuffer, command);
    }

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, _nodeDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, _nodeDataBuffer);

    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, _drawCommandBuffer);

    _shader.bind();
    _mesh.bind();

    const glm::mat4 viewProj = camera.getProjectionMatrix() * camera.getViewMatrix();
    _shader.set("modelView", camera.getViewMatrix());
    _shader.set("modelViewProj", viewProj);

    glDisable(GL_DEPTH_TEST);
    glDrawElementsIndirect(GL_TRIANGLES, GL_UNSIGNED_INT, nullptr);
    glEnable(GL_DEPTH_TEST);
}
}