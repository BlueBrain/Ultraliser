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
#include <data/structures/List.h>
#include <algorithms/utilities/Sorting.h>

namespace Ultraliser
{

// Create a new node containing 'd', and link it to       //
// 'p' on the left (prev) and to 'n' on the right (next). //

Node::Node(const Node *p, const void *d, const Node *n)
{
    data=(void *)d;
    if ((_prevNode=(Node *)p) != nullptr) _prevNode->_nextNode = this;
    if ((_nextNode=(Node *)n) != nullptr) _nextNode->_prevNode = this;
}


// Destroy and unlink the node

Node::~Node()
{
    if (_prevNode != nullptr) _prevNode->_nextNode = _nextNode;
    if (_nextNode != nullptr) _nextNode->_prevNode = _prevNode;
}


/////////// Constructor from list ///////////////////

List::List(const void **d, int n)
{
    _lHead = _lTail = nullptr; _lNumElements = 0;
    for (int i=0; i<n; i++) appendTail(d[i]);
}

///////////////////////// Destructor //////////////////////////

List::~List()
{
    while (_lHead != nullptr) removeCell(_lHead);
}

////////////////// Append an element //////////////////

void List::appendHead(const void *d)
{
    _lHead = new Node(nullptr, d, _lHead);
    if (_lTail == nullptr) _lTail = _lHead;
    _lNumElements++;
}

void List::appendTail(const void *d)
{
    _lTail = new Node(_lTail, d, nullptr);
    if (_lHead == nullptr) _lHead = _lTail;
    _lNumElements++;
}

void List::insertAfter(Node *b, const void *d)
{
    Node *nn = new Node(b, d, b->next());
    if (b == _lTail) _lTail = nn;
    _lNumElements++;
}

////////////////// Appends a list //////////////////

void List::appendList(const List *l)
{
    Node *n = l->_lTail;

    while (n != nullptr)
    {
        appendHead(n->data);
        n=n->prev();
    }
}

////////////////// Joins a list to the l_tail //////////////////

void List::joinTailList(List *l)
{
    if (l->_lNumElements == 0) return;
    if (_lTail != nullptr)
    {
        _lTail->_nextNode = l->_lHead; l->_lHead->_prevNode = _lTail; _lTail = l->_lTail;
        _lNumElements += l->_lNumElements;
    }
    else
    {
        _lHead = l->_lHead; _lTail = l->_lTail; _lNumElements = l->_lNumElements;
    }
    l->_lHead = l->_lTail = nullptr; l->_lNumElements = 0;
}

//// Moves node 'n' from this list to the end of 'l'. \n O(1).
void List::moveNodeTo(Node *n, List *l)
{
    Node *pn = n->_prevNode, *nn = n->_nextNode;
    n->_prevNode = l->_lTail; n->_nextNode = nullptr;

    if (l->_lNumElements) l->_lTail->_nextNode = n; else l->_lHead = n;
    l->_lTail = n;
    l->_lNumElements++;

    _lNumElements--;
    if (pn != nullptr) pn->_nextNode = nn; else _lHead = nn;
    if (nn != nullptr) nn->_prevNode = pn; else _lTail = pn;
}

//// Removes the first node and returns the corresponding data /////

void *List::popHead()
{
    void *data = (_lHead != nullptr)?(_lHead->data):(nullptr);
    if (_lHead != nullptr) removeCell(_lHead);
    return data;
}

//// Removes the last node and returns the corresponding data /////

void *List::popTail()
{
    void *data = (_lTail != nullptr)?(_lTail->data):(nullptr);
    if (_lTail != nullptr) removeCell(_lTail);
    return data;
}

//////////////////// Removes an element //////////////////

int List::removeNode(const void *d)
{
    Node *tmp = _lHead;
    int i=1;

    while (tmp != nullptr)
        if (tmp->data == d)
        {
            removeCell(tmp);
            return i;
        }
        else {tmp=tmp->_nextNode; i++;}

    return 0;
}


//////////////////// Removes an element //////////////////

int List::removeNode(int i)
{
    Node *tmp = _lHead;

    while (tmp!=nullptr && i--) tmp=tmp->_nextNode;
    if (tmp==nullptr) return 0;

    removeCell(tmp);
    return 1;
}


//////////////////// Gets a node //////////////////

Node *List::getNode(int i) const
{
    Node *tmp = _lHead;

    while (tmp!=nullptr && i--) tmp=tmp->_nextNode;
    return tmp;
}


//////////////////// Removes a node //////////////////

void List::removeCell(Node *n)
{
    if (n==_lHead) _lHead = n->_nextNode;
    if (n==_lTail) _lTail = n->_prevNode;
    delete(n);
    _lNumElements--;
}


////////////////// Garbage collection //////////////

void List::freeCell(Node *n)
{
    free(n->data);
    removeCell(n);
}

void List::freeNode(void *d)
{
    free(d);
    removeNode(d);
}

//////////////////// Belonging check /////////////////

Node *List::containsNode(const void *d) const
{
    Node *tmp = _lHead;

    while (tmp != nullptr)
        if (tmp->data == d) return tmp;
        else tmp=tmp->_nextNode;

    return nullptr;
}

//////////////////// Replaces a node /////////////////

Node *List::replaceNode(const void *od, const void *nd)
{
    Node *tmp = containsNode(od);
    if (tmp != nullptr) {tmp->data = (void *)nd; return tmp;}
    appendTail(nd);
    return _lTail;
}

//////////////////////// Garbage collector /////////////////////

void List::freeNodes()
{
    while (_lHead != nullptr) freeCell(_lHead);
}

//////////////////////// Garbage collector /////////////////////

void List::removeNodes()
{
    while (_lHead != nullptr) removeCell(_lHead);
}


///// Conversion to array ///////

void **List::toArray() const
{
    Node *n = _lHead;
    int i;
    void **array;

    if (_lNumElements == 0) return nullptr;
    array = (void **)malloc(sizeof(void *)*_lNumElements);
    if (array == nullptr) return nullptr;
    for (i=0; i<_lNumElements; i++, n=n->_nextNode) array[i] = n->data;

    return array;
}

///// Sorts the list /////////

int List::sort(int (*comp)(const void *, const void *))
{
    void **array;
    int ne = _lNumElements-1;

    if (_lNumElements < 2) return 0;
    if ((array = toArray()) == nullptr) return 1;

    jqSort(array, _lNumElements, comp);
    removeNodes();
    for (; ne >= 0; ne--) appendHead(array[ne]);
    free(array);

    return 0;
}

}
