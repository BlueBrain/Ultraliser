/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marco Attene < IMATI-GE / CNR >
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
 *
 * The content of this file is based on MeshFix. The code has been modified under the terms of
 * the GNU General Public License as published by the Free Software Foundation either version 3 of
 * the License, or (at your option) any later version.
 * MeshFix has a dual license for free and commercial use. For further information, please refer
 * to the original repository at < https://github.com/MarcoAttene/MeshFix-V2.1>.
 **************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "Heap.h"

namespace Ultraliser
{

AbstractHeap::AbstractHeap(int size)
{
    _heap = new void *[size + 1];
    _numberElements = 0;
    _maxNumberElements = size;
    _positions = nullptr;
}

AbstractHeap::~AbstractHeap()
{
    delete(_heap);
}

int AbstractHeap::_upHeap(int k)
{
    if (k < 2) return k;

    void *t = _heap[k];
    int fk = (k%2)?((k-1)/2):(k/2);
    void *f = _heap[fk];

    if (_compare(t, f) <= 0)
    {
        _heap[k] = f;
        _heap[fk] = t;
        if (_positions != nullptr)
        {
            int* ff = (int*) f;
            int* tt = (int*) t;

            _positions[*ff] = k;
            _positions[*tt] = fk;
        }
        return _upHeap(fk);
    }
    return k;
}

int AbstractHeap::_downHeap(int k)
{
    int j;

    void *t = _heap[k];
    int fk = (_numberElements%2)?((_numberElements-1)/2):(_numberElements/2);
    if (k > fk) return k;

    j = k+k;
    if (j < _numberElements && _compare(_heap[j], _heap[j+1]) >= 0) j++;
    void *f = _heap[j];
    if (_compare(t, f) >= 0)
    {
        _heap[k] = f;
        _heap[j] = t;
        if (_positions != nullptr)
        {
            int* ff = (int*) f;
            int* tt = (int*) t;
            _positions[*ff] = k;
            _positions[*tt] = j;
        }
        return _downHeap(j);
    }

    return k;
}

int AbstractHeap::insert(void *t)
{
    if (_numberElements == _maxNumberElements) return -1;

    _heap[++_numberElements] = t;


    int* tt = (int*) t;
    if (_positions != nullptr) _positions[*tt] = _numberElements;
    return _upHeap(_numberElements);
}

void *AbstractHeap::removeHead()
{
    void *t = _heap[1];
    int* tt = (int*) t;
    if (_positions != nullptr) _positions[*tt] = 0;
    _heap[1] = _heap[_numberElements--];
    if (_numberElements)
    {
        int* hh = (int*) _heap[1];
        if (_positions != nullptr) _positions[*hh] = 1;
        _downHeap(1);
    }

    return t;
}

} 
