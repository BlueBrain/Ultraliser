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

#include "VolumeMaterial.h"

#include "../Model.h"

namespace scr
{
  VolumeMaterial::VolumeMaterial(Shader* shader, Texture3D* texture)
   : Material(shader)
   , _volume(texture)
  {
    if(_volume == nullptr)
      throw std::runtime_error("VolumeMaterial: Null 3D Texture passed");

    _volume->makeHandleResident();
  }

  void VolumeMaterial::setBounds(const glm::vec3& min, const glm::vec3& max)
  {
    _minBound = min;
    _maxBound = max;
  }

  void VolumeMaterial::render(const Camera& camera, const Model& model)
  {
    (void)model;

    glEnable(GL_DEPTH_TEST);

    const glm::mat4 viewProj = camera.getProjectionMatrix() * camera.getViewMatrix();
    const glm::mat4 invView = glm::inverse(camera.getViewMatrix());
    _shader->set("modelView", camera.getViewMatrix());
    _shader->set("modelViewProj", viewProj);
    _shader->set("invModelView", invView);

    const glm::vec3 lightDir (-1.f, 1.f, 1.f);
    _shader->set("lightDirection", lightDir);

    _shader->set("localMinBound", _minBound);
    _shader->set("localMaxBound", _maxBound);

    _shader->set("volumeTexture", _volume->getBindlessHandle());

    glDrawElements(GL_TRIANGLES, 12 * 3, GL_UNSIGNED_INT, (void*)0);
  }
}