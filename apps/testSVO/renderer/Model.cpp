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

#include "Model.h"

#include <stdexcept>

namespace scr
{
  Model::Model(CubeMesh* mesh,
               Material* material)
   : _mesh(mesh)
   , _material(material)
  {
    if(_mesh == nullptr)
      throw std::runtime_error("Model: Null mesh passed");

    if(_material == nullptr)
      throw std::runtime_error("Model: Null material passed");
  }

  Model::~Model()
  {
    unbind();
    _mesh->destroy();
    _material->destroy();
  }

  void Model::bind() const
  {
    _mesh->bind();
    _material->bind();
  }

  void Model::unbind() const
  {
    _mesh->unbind();
    _material->unbind();
  }

  void Model::renderModel(const Camera& camera)
  {
    _material->render(camera, *this);
  }
}