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

#ifndef SIMCRUSHERRENDERER_CAMERA_H
#define SIMCRUSHERRENDERER_CAMERA_H

#include <glm/glm.hpp>

namespace svorender
{
class Camera
{
public:
    Camera(float nearPlane, float farPlane, float fov);

    void setPosition(const glm::vec3 &position);
    void updatePosition(const glm::vec3 &deltaPos);
    void setRotation(const glm::vec3 &axisAngles);
    void updateRotation(const glm::vec3 &deltaAxisAngles);

    const glm::mat4 &getViewMatrix() const;
    const glm::mat4 &getProjectionMatrix() const;

    void onScreenResize(uint32_t width, uint32_t height);
    void onZoom(float val);

    void update();

private:
    void updateProjectionMatrix();

private:
    float _near;
    float _far;
    float _fov;

    float _aspectRatio;

    glm::vec3 _position;
    glm::vec3 _rotation;

    glm::mat4 _viewMatrix;
    glm::mat4 _projectionMatrix;

    bool _dirty;
};
}

#endif