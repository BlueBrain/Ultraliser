/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#ifndef ULTRALISER_MESHES_SIMPLE_PRIMITIVES_TEXTURE_HH
#define ULTRALISER_MESHES_SIMPLE_PRIMITIVES_TEXTURE_HH

#include <math/Math.h>
#include <data/meshes/simple/primitives/Quad.hh>

namespace Ultraliser
{

/**
 * @brief QuadPtr
 */
typedef SimpleQuad* QuadPtr;

/**
 * @brief Quads
 */
typedef std::vector< SimpleQuad > Quads;

/**
 * @brief Texture
 */
typedef Vector2f SimpleTexture;

/**
 * @brief TexturePtr
 */
typedef SimpleTexture* TexturePtr;

/**
 * @brief Textures
 */
typedef std::vector< SimpleTexture > Textures;

/**
 * @brief TextureID
 */
typedef Vec3i_64 TextureID;

/**
 * @brief TextureIDs
 */
typedef std::vector <TextureID> TextureIDs;

}

#endif // ULTRALISER_MESHES_SIMPLE_PRIMITIVES_TEXTURE_HH