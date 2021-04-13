/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#include "Section.h"

namespace Ultraliser
{

Section::Section(const uint64_t &index)
    : _index(index)
{
    /// EMPTY CONSTRUCTOR
}

uint64_t Section::getIndex() const
{
    return _index;
}

void Section::addParentIndex(const uint64_t index)
{
    _parentsIndices.push_back(index);
}

void Section::addChildIndex(const uint64_t index)
{
    _childrenIndices.push_back(index);
}

void Section::addSample(Sample* sample)
{
    _samples.push_back(sample);
}

Samples Section::getSamples() const
{
    return _samples;
}

std::vector< uint64_t > Section::getParentIndices() const
{
    return _parentsIndices;
}

std::vector< uint64_t > Section::getChildrenIndices() const
{
    return _childrenIndices;
}

}
