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

#include "VolumeTreeModel.h"

namespace
{
class ShaderFactory
{
public:
    static svorender::Shader create()
    {
        svorender::ShaderConfig volumeShader;
        volumeShader.vertexShader = "./shaders/volrenderV.glsl";
        volumeShader.fragmentShader = "./shaders/volrenderF.glsl";
        return svorender::Shader(volumeShader);
    }
};
}

namespace svorender
{
VolumeTreeModel::VolumeTreeModel(const glm::vec3 &min, const glm::vec3 &max, Texture3D volume)
    : _min(min)
    , _max(max)
    , _volume(std::move(volume))
    , _mesh(min, max)
    , _shader(ShaderFactory::create())
{
    _volume.makeHandleResident();
}

void VolumeTreeModel::render(const Camera &camera)
{
    _mesh.bind();
    _shader.bind();

    auto viewProj = camera.getProjectionMatrix() * camera.getViewMatrix();
    auto invView = glm::inverse(camera.getViewMatrix());
    _shader.set("modelView", camera.getViewMatrix());
    _shader.set("modelViewProj", viewProj);
    _shader.set("invModelView", invView);

    auto lightDir = glm::normalize(glm::vec3(-1.f, 1.f, 1.f));
    _shader.set("lightDirection", lightDir);

    _shader.set("localMinBound", _min);
    _shader.set("localMaxBound", _max);

    _shader.set("volumeTexture", _volume.getBindlessHandle());

    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, (void *)0);
}
}