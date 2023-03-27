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

#pragma once

#include "../CubeMesh.h"
#include "../Model.h"
#include "../Shader.h"

#include <vector>

namespace svorender
{
class NodeTreeModel : public Model
{
public:
    NodeTreeModel();

    void addNode(const glm::vec3 &position, float scale);
    void render(const Camera &camera) override;

private:
    CubeMesh _mesh;
    Shader _shader;
    GLuint _nodeDataBuffer = GL_INVALID_VALUE;
    GLuint _drawCommandBuffer = GL_INVALID_VALUE;

    std::vector<glm::vec4> _nodes;
    bool _dirty = false;
};
}