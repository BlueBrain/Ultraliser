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

#include "BoundRenderMaterial.h"

#include "../Model.h"

namespace scr
{
  BoundRenderMaterial::BoundRenderMaterial(Shader* shader)
   : Material(shader)
  {
  }

  void BoundRenderMaterial::render(const Camera& camera, const Model& model)
  {
    (void)model;
    glDisable(GL_DEPTH_TEST);
    const glm::mat4 viewProj = camera.getProjectionMatrix() * camera.getViewMatrix();
    _shader->set("modelView", camera.getViewMatrix());
    _shader->set("modelViewProj", viewProj);

    glDrawElements(GL_TRIANGLES, 12 * 3, GL_UNSIGNED_INT, (void*)0);
  }
}