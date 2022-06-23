/***************************************************************************************************
 * Copyright (c) 2016 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#pragma once

#include "Volume.h"
#include <utilities/Timer.h>
#include <data/volumes/grids/Grids.h>

template
void writeNRRD< uint8_t >(const std::string &, const VolumeGridU8*);

template
void writeNRRD< uint16_t >(const std::string &, const VolumeGridU16*);

template
void writeNRRD< uint32_t >(const std::string &, const VolumeGridU32*);

template
void writeNRRD< uint64_t >(const std::string &, const VolumeGridU64*);

template
void writeNRRD< float >(const std::string&, const VolumeGridF32*);

template
void writeNRRD< double >(const std::string&, const VolumeGridF64*);

template
void writeRAW< uint8_t >(const std::string&, const VolumeGridU8*);

template
void writeRAW< uint16_t >(const std::string&, const VolumeGridU16*);

template
void writeRAW< uint32_t >(const std::string&, const VolumeGridU32*);

template
void writeRAW< uint64_t >(const std::string&, const VolumeGridU64*);

template
void writeVOL< uint8_t >(const std::string &, const VolumeGridU8*);

template
void writeVOL< uint16_t >(const std::string &, const VolumeGridU16*);

template
void writeVOL< uint32_t >(const std::string &, const VolumeGridU32*);

template
void writeVOL< uint64_t >(const std::string &, const VolumeGridU64*);

template
void writeVOL< float >(const std::string&, const VolumeGridF32*);

template
void writeVOL< double >(const std::string&, const VolumeGridF64*);
