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

#include "Texture3D.h"
#include "Common.h"

#include <stdexcept>

namespace scr
{
  Texture3D::Texture3D(const uint32_t dimX, 
                       const uint32_t dimY,
                       const uint32_t dimZ,
                       const std::vector<uint8_t>& volumeData)
   : _gridX(dimX)
   , _gridY(dimY)
   , _gridZ(dimZ)
   , _oglTexture(GL_INVALID_VALUE)
  {
    if(volumeData.size() < _gridX * _gridY * _gridZ)
      throw std::runtime_error("Texture3D: Given volume data is smaller than grid measures");

    glGenTextures(1, &_oglTexture);
    if(_oglTexture == GL_INVALID_VALUE)
      throw std::runtime_error("Texture3D: Could not create texture");

    bind();
    
    glTexStorage3D(GL_TEXTURE_3D, 1, GL_R8, _gridX, _gridY, _gridZ);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, _gridX, _gridY, _gridZ, GL_RED, GL_UNSIGNED_BYTE, volumeData.data());

    //glGenerateMipmap(GL_TEXTURE_3D);

    // Generate bindless handle
    _oglTextureHandle = glGetTextureHandleARB(_oglTexture);

    if(_oglTextureHandle == GL_INVALID_VALUE)
      throw std::runtime_error("Texture3D: Could not get bindless handle for texture");
  }

  void Texture3D::bind() const
  {
    glBindTexture(GL_TEXTURE_3D, _oglTexture);
  }

  void Texture3D::unbind() const
  {
    glBindTexture(GL_TEXTURE_3D, 0);
  }

  void Texture3D::makeHandleResident() const
  {
    glMakeTextureHandleResidentARB(_oglTextureHandle);
  }

  void Texture3D::makeHandleNonResident() const
  {
    glMakeTextureHandleNonResidentARB(_oglTextureHandle);
  }

  GLuint64 Texture3D::getBindlessHandle() const
  {
    return _oglTextureHandle;
  }

  void Texture3D::destroy()
  {
    makeHandleNonResident();
    glDeleteTextures(1, &_oglTexture);
  }
/*
  std::vector<uint8_t> Texture3D::octreeToVolume(const sc::SparseOctree& tree) const
  {
    const glm::vec3 size = _maxBound - _minBound;

    const sc::Point3DF treeMinBound = tree.getBounds().min;

    // Cell size
    const float dx = size.x / static_cast<float>(_gridResolution);
    const float dy = size.y / static_cast<float>(_gridResolution);
    const float dz = size.z / static_cast<float>(_gridResolution);

    // Cell center
    const float centerx = dx / 2.f;
    const float centery = dy / 2.f;
    const float centerz = dz / 2.f;

    const uint32_t frameLen = _gridResolution * _gridResolution;

    const sc::VoxelPointTest algorithm;

    std::vector<uint8_t> result;
    result.resize(_gridResolution 
                  * _gridResolution 
                  * _gridResolution, 
                  0u);
    
    for(size_t i = 0; i < _gridResolution; i++)
    {
      for(size_t j = 0; j < _gridResolution; j++)
      {
        for(size_t k = 0; k < _gridResolution; k++)
        {
          const sc::Point3DF samplePoint = treeMinBound 
                                            + sc::Point3DF(dx * k + centerx,
                                                           dy * j + centery,
                                                           dz * i + centerz);

          const uint8_t sample = sc::sampleShape(tree, 
                                                 algorithm, 
                                                 sc::PointShape(samplePoint));

          result[frameLen * i + _gridResolution * j + k] = sample;
        }
      }
    }

    return result;
  }
  */
}