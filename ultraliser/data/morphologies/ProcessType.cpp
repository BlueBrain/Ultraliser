/***************************************************************************************************
 * Copyright (c) 2016 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#include "ProcessType.h"

int64_t mapNeuronProcessTypeToSWCIndex(const PROCESS_TYPE& type)
{
    switch (type)
    {
    case PROCESS_TYPE::SOMA:
        return 1;
        break;

    case PROCESS_TYPE::NEURON_AXON:
        return 2;
        break;

    case PROCESS_TYPE::NEURON_BASAL_DENDRITE:
        return 3;
        break;

    case PROCESS_TYPE::NEURON_APICAL_DENDRITE:
        return 4;
        break;

    default:
        return 0;
        break;
    }
}

PROCESS_TYPE mapNeuronH5IndexToType(const size_t& index)
{
    if (index == 0)
        return PROCESS_TYPE::UNKNOWN_PROCESS;
    else if (index == 1)
        return PROCESS_TYPE::SOMA;
    else if (index == 2)
        return PROCESS_TYPE::NEURON_AXON;
    else if (index == 3)
        return PROCESS_TYPE::NEURON_BASAL_DENDRITE;
    else if (index == 4)
        return PROCESS_TYPE::NEURON_APICAL_DENDRITE;
    else
        return PROCESS_TYPE::UNKNOWN_PROCESS;
}
