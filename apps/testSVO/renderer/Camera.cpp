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

#include "Camera.h"

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace svorender
{
Camera::Camera(float nearPlane, float farPlane, float fov)
    : _near(nearPlane)
    , _far(farPlane)
    , _fov(fov)
    , _aspectRatio(1.f)
    , _position(glm::vec3(0.f))
    , _rotation(glm::vec3(0.f))
    , _viewMatrix(glm::mat4(1.f))
    , _projectionMatrix(glm::mat4(0.f))
    , _dirty(true)
{
    updateProjectionMatrix();
}

void Camera::setPosition(const glm::vec3 &position)
{
    _position = position;
    _dirty = true;
}

void Camera::updatePosition(const glm::vec3 &deltaPos)
{
    _position += deltaPos;
    _dirty = true;
}

void Camera::setRotation(const glm::vec3 &axisAngles)
{
    _rotation = axisAngles;
    _dirty = true;
}

void Camera::updateRotation(const glm::vec3 &deltaAxisAngles)
{
    _rotation += deltaAxisAngles;
    _dirty = true;
}

const glm::mat4 &Camera::getViewMatrix() const
{
    return _viewMatrix;
}

const glm::mat4 &Camera::getProjectionMatrix() const
{
    return _projectionMatrix;
}

void Camera::onScreenResize(uint32_t width, uint32_t height)
{
    auto fWidth = static_cast<float>(width);
    auto fHeight = static_cast<float>(height);
    _aspectRatio = fWidth / fHeight;
    updateProjectionMatrix();
}

void Camera::onZoom(float val)
{
    _fov += val;
    _fov = std::min(45.0f, std::max(1.0f, _fov));
    updateProjectionMatrix();
}

void Camera::update()
{
    if (!_dirty)
    {
        return;
    }

    _dirty = false;

    glm::quat yawQ = glm::quat(glm::vec3(0.0f, _rotation.y, 0.0f));
    yawQ = glm::normalize(yawQ);
    glm::mat4 yaw = glm::mat4_cast(yawQ);

    glm::quat pitchQ = glm::quat(glm::vec3(_rotation.x, 0.0f, 0.0f));
    pitchQ = glm::normalize(pitchQ);
    glm::mat4 pitch = glm::mat4_cast(pitchQ);

    _viewMatrix = glm::translate(glm::mat4(1.f), _position) * yaw * pitch;
}

void Camera::updateProjectionMatrix()
{
    auto rads = glm::radians(_fov);
    auto tanRads = glm::tan(rads);
    _projectionMatrix = glm::mat4(0.f);
    _projectionMatrix[0].x = 1.f / (_aspectRatio * tanRads);
    _projectionMatrix[1].y = 1.f / tanRads;
    _projectionMatrix[2].z = (_far + _near) / (_near - _far);
    _projectionMatrix[3].z = 2.f * _near * _far / (_near - _far);
    _projectionMatrix[2].w = -1.f;
}
}