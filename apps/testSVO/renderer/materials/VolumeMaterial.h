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

#ifndef SIMCRUSHERRENDERER_VOLUMEMATERIAL_H
#define SIMCRUSHERRENDERER_VOLUMEMATERIAL_H

#include "../Material.h"
#include "../Texture3D.h"

#include <stdexcept>

namespace scr
{
class VolumeMaterial : public Material
{
public:
    VolumeMaterial(Shader *shader, Texture3D *texture);
    void setBounds(const glm::vec3 &min, const glm::vec3 &max);

    void render(const Camera &camera, const Model &model);

private:
    Texture3D *_volume;
    glm::vec3 _minBound;
    glm::vec3 _maxBound;
};
}

#endif