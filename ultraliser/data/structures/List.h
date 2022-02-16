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

#ifndef ULTRALISER_DATA_STRUCTURES_LIST_H
#define ULTRALISER_DATA_STRUCTURES_LIST_H

#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The Node class
 * Generic node of a doubly likned list.
 */
class Node
{
    /// This is to make methods in 'List' able to modify _prevNode and _nextNode
    friend class List;

public:

    /**
     * @brief data
     * Actual data stored in the node.
     */
    void* data;

    /**
     * @brief Node
     * Creates an isolated node storing 'd'.
     * @param d
     */
    Node(const void *d)
    {
        data = (void*) d;
        _prevNode = _nextNode = nullptr;
    }

    /**
     * @brief Node
     * Creates a new node storing 'd' and links it to a previous node 'p'
     * and to a next one 'n'.
     * @param p
     * @param d
     * @param n
     */
    Node(const Node *p, const void *d, const Node *n);
    ~Node();

    /**
     * @brief prev
     * Returns the previous node in the list, possibly nullptr.
     * @return
     */
    inline Node *prev() const
    {
        return _prevNode;
    }

    /**
     * @brief next
     * Returns the next node in the list, possibly nullptr
     * @return
     */
    inline Node *next() const
    {
        return _nextNode;
    }

protected:

    /**
     * @brief _prevNode
     *
     */
    Node* _prevNode;

    /**
     * @brief _nextNode
     */
    Node* _nextNode;
};

/**
 * @brief The List class
 * Doubly linked list.
 */
class List
{ 
protected :

    /**
     * @brief _lHead
     * First node pointer
     */
    Node *_lHead;

    /**
     * @brief _lTail
     * Last node pointer
     */
    Node *_lTail;

    /**
     * @brief _lNumElements
     * Number of elements in the list
     */
    uint64_t _lNumElements;

public :

    /**
     * @brief List
     * Creates an empty list
     */
    List()
    {
        _lHead = _lTail = nullptr;
        _lNumElements = 0;
    }

    /**
     * @brief List
     * Creates a list containing an element 'd' (singleton)
     * @param d
     */
    List(const void *d)
    {
        _lHead = _lTail = new Node(d);
        _lNumElements = 1;
    }

    /**
     * @brief List
     * Creates a list out of an array 'd' made of 'n' elements.
     * @param d
     * @param n
     */
    List(const void **d, int n);

    /**
     * @brief List
     * Creates a duplicated list.
     * @param l
     */
    List(List& l)
    {
        _lHead = _lTail = nullptr;
        _lNumElements = 0;
        appendList(&l);
    }

    /**
     * @brief List
     * Creates a duplicated list.
     * @param l
     */
    List(List* l)
    {
        _lHead = _lTail = nullptr;
        _lNumElements = 0;
        appendList(l);
    }
    ~List();

    /**
     * @brief head
     * Gets the first node, nullptr if empty. O(1).
     * @return
     */
    Node *head() const { return _lHead; }

    /**
     * @brief tail
     * Gets the last node, nullptr if empty. O(1).
     * @return
     */
    Node *tail() const { return _lTail; }

    /**
     * @brief numberElements
     * Gets the number of elements. O(1).
     * @return
     */
    uint64_t numberElements() const { return _lNumElements; }

    /**
     * @brief appendHead
     * Appends a new node storing 'd' to the head.
     * Complexity O(1).
     * @param d
     */
    void appendHead(const void *d);

    /**
     * @brief appendTail
     * Appends a new node storing 'd' to the tail.
     * Complexity O(1).
     * @param d
     */
    void appendTail(const void *d);

    /**
     * @brief insertAfter
     * Inserts a new node storing 'd' right after 'n'.
     * Complexity O(1).
     * @param n
     * @param d
     */
    void insertAfter(Node *n, const void *d);

    /**
     * @brief removeNode
     *  Deletes and removes the node containing 'd'.
     * Returns its position, 0 if 'd' is not in the list.
     * Complexity O(N).
     * @param d
     * @return
     */
    int removeNode(const void *d);

    /**
     * @brief removeNode
     * Deletes and i'th node (starting from 0).
     * Returns 0 if the list has less than i + 1 nodes.
     * Complexity O(N).
     * @param i
     * @return
     */
    int removeNode(int i);

    /**
     * @brief getNode
     * Returns the node at position 'i' (starting from 0).
     * Returns nullptr if the list has less than i + 1 nodes.
     * Complexity O(N).
     * @param i
     * @return
     */
    Node *getNode(int i) const;

    /**
     * @brief removeCell
     * Deletes and removes the node 'n' from the list.
     * Complexity O(1).
     * @param n
     */
    void removeCell(Node *n);

    /**
     * @brief appendList
     * Appends a list 'l' to the head by duplicating nodes in 'l'.
     * Complexity O(l->Number Elements).
     * @param l
     */
    void appendList(const List *l);

    /**
     * @brief joinTailList
     * Appends a list 'l' to the tail by linking the first node of 'l' to the
     * last one of this list. 'l' becomes empty.
     * Complexity O(1).
     * @param l
     */
    void joinTailList(List *l);

    /**
     * @brief moveNodeTo
     *  Moves node 'n' from this list to the end of 'l'.
     * Complexity O(1).
     * @param n
     * @param l
     */
    void moveNodeTo(Node *n, List *l);

    /**
     * @brief popHead
     * Deletes and removes the first node. Returns its data.
     * Complexity O(1).
     * @return
     */
    void *popHead();

    /**
     * @brief popTail
     * Deletes and removes the last node. Returns its data.
     * Complexity O(1).
     * @return
     */
    void *popTail();

    /**
     * @brief freeCell
     * Deletes and removes the node 'n' from the list and frees data memory.
     * Complexity O(1).
     * Warning. This method uses the free() function to to dispose the memory
     * space used by the data stored in the node. This means that such data
     * should have been allocated through malloc(), calloc() or realloc(),
     * and not through the 'new' operator.
     * On some systems, however, the 'delete' operator simply calls 'free()'
     * right after the execution of the proper object destructor so, if the
     * object does not need to free internally allocated memory, it is safe to
     * dispose the memory trhough free() although the object was allocated by
     * 'new'. This works on Linux Fedora Core 2 distributions.
     * @param n
     */
    void freeCell(Node *n);

    /**
     * @brief freeNode
     * Deletes and removes the node storing 'd' and frees the memory occupied
     * by 'd' itself. O(N).
     * Warning. Read the comment for the method 'freeCell()'
     * @param d
     */
    void freeNode(void *d);

    /**
     * @brief containsNode
     * Returns the node storing 'd'. nullptr if not found.
     * Complexity O(N).
     * @param d
     * @return
     */
    Node *containsNode(const void *d) const;

    /**
     * @brief replaceNode
     * Replaces old_n with new_n. The Node containing new_n is returned.
     * Complexity O(N).
     * @param old_n
     * @param new_n
     * @return
     */
    Node *replaceNode(const void *old_n, const void *new_n);

    /**
     * @brief freeNodes
     * Deletes and removes all the nodes and frees data memory.
     * Complexity O(N).
     * Warning. Read the comment for the method 'freeCell()'
     */
    void freeNodes();

    /**
     * @brief removeNodes
     * Deletes and removes all the nodes. O(N).
     */
    void removeNodes();

    /**
     * @brief toArray
     * Creates an array out of the list. O(N).
     * @return
     */
    void **toArray() const;

    /**
     * @brief sort
     * Sorts the list using 'comp' as comparison function for two elements.
     * Complexity O(N^2).
     * This method uses the QuickSort algorithm for sorting, thus the complexity
     * is N^2 in the worst case, but it is actually much faster in the
     * average case. If, however, there is the need to have a guaranteed
     * O(NlogN) complexity, it is possible to implement a heap based on the
     * 'abstractHeap' class. See the documentation of the standard 'qsort' library
     * function for details on the prototype of the comparison function 'comp'.
     * @return
     */
    int sort(int (*comp)(const void *, const void *));
};

/// Convenience macro to scan the nodes of a list.
#define FOR_EACH_NODE(l, n)                                                     \
    for ((n) = (l).head(); (n) != nullptr; (n) = (n)->next())

/// Convenience macro to circulate around the nodes of a list 'l' starting from
/// node 'm'. Must exit with break or return.
#define FOR_EACH_NODE_CIRCULAR(l, m, n)                                         \
    for ((n) = (m); ; (n) = ((n) != (l).tail())?((n)->next()):((l).head()))
}
#endif // ULTRALISER_DATA_STRUCTURES_LIST_H

