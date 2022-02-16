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

#ifndef ULTRALISER_DATA_STRUCTURES_HEAP_H
#define ULTRALISER_DATA_STRUCTURES_HEAP_H

#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The AbstractHeap class
 * AbstractHeap is the base class for implementing heaps.
 *
 * Each implementation (class extension) must define the method
 * compare to be used for sorting the heap. If the objects being
 * sorted are non-negative numbers, a special implementation may
 * use the field positions to record the index of each element
 * within the heap. This feature is useful when there is a need to
 * re-sort an element whose cost changes after its insertion into
 * the heap. The array 'positions' must be allocated by
 * the extended class constructor, and must be able to contain NMAX+1
 * integer numbers, where NMAX is the maximum value that can be
 * assumed by an object.
 */
class AbstractHeap
{
protected:

    /**
     * @brief _heap
     * Heap data is stored here
     */
    void** _heap;

    /**
     * @brief _numberElements
     * Current number of elements
     */
    int _numberElements;

    /**
     * @brief _maxNumberElements
     * Maximum number of elements
     */
    int _maxNumberElements;

    /**
     * @brief _positions
     * Optional pointer to an array of positions
     */
    int *_positions;

    /**
     * @brief _upheap
     * Moves the i'th object up on the heap
     * @param i
     * @return
     */
    int _upHeap(int i);

    /**
     * @brief _downHeap
     * Moves the i'th object down on the heap
     * @param i
     * @return
     */
    int _downHeap(int i);

    /**
     * @brief _compare
     * Comparison of two heap elements.
     * This function must be implemented in the extended class.
     * The return value must be <0 if a<b, >0 if a>b or 0 if a=b.
     * @param a
     * @param b
     * @return
     */
    virtual int _compare(const void *a, const void *b) = 0;

public :

    /**
     * @brief AbstractHeap
     * Creates a heap which can contain up to 'n' elements.
     * @param n
     */
    AbstractHeap(int n);

    /**
     * @brief ~AbstractHeap
     * Default destructor.
     */
    virtual ~AbstractHeap() = 0;

    /**
     * @brief insert
     * Inserts 'e' into the heap.
     * Inserts an element 'e' into the heap in the correct position, according
     * to the method compare. If the insertion fails because the heap is
     * full, -1 is returned, otherwise the index position of the newly inserted
     * element is returned.
     * @param e
     * @return
     */
    int insert(void *e);

    /**
     * @brief isEmpty
     * Returns TRUE if the heap is empty.
     * @return
     */
    int isEmpty() const {return (_numberElements == 0);}

    /**
     * @brief getHead
     * Returns the first element of the heap.
     * @return
     */
    void *getHead() const {return _heap[1];}

    /**
     * @brief removeHead
     * Removes and returns the first element after rearranging the heap.

     * @return
     */
    void *removeHead();

    /**
     * @brief flush
     * Removes all the elements.
     */
    void flush() { _numberElements=0; }
};

}

#endif // ULTRALISER_DATA_STRUCTURES_HEAP_H
