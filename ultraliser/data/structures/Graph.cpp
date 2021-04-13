/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#include "Graph.h"

namespace Ultraliser
{

GraphEdge *GraphNode::getEdge(GraphNode *gn)
{
    GraphEdge *ge;
    Node *n = edges.head();
    while (n!=nullptr)
        if ((ge=((GraphEdge *)n->data))->oppositeNode(this)==gn) return ge;
        else n=n->next();
    return nullptr;
}

GraphEdge::GraphEdge(GraphNode *a, GraphNode *b)
{
    n1=a; n2=b;
    n1->edges.appendHead(this);
    n2->edges.appendHead(this);
}


GraphEdge *Graph::createEdge(GraphNode *n1, GraphNode *n2)
{
    Node *n;
    FOR_EACH_NODE(n1->edges, n)
            if (((GraphEdge *)n->data)->hasNode(n2))
            return (GraphEdge *)n->data;

    edges.appendHead(new GraphEdge(n1, n2));
    return (GraphEdge *)edges.head()->data;
}


void GraphEdge::collapse()
{
    Node *n;
    GraphEdge *e;
    GraphNode *nx;

    while ((e = (GraphEdge *)n2->edges.popHead()) != nullptr)
        if (e != this)
        {
            ((e->n1 == n2)?(e->n1):(e->n2)) = n1;
            n1->edges.appendHead(e);
        }

    FOR_EACH_NODE(n1->edges, n)
    {
        e = (GraphEdge *)n->data;
        if (!e->isUnlinked()) e->oppositeNode(n1)->mask = 0;
    }

    n2->mask = 1;
    FOR_EACH_NODE(n1->edges, n)
    {
        e = (GraphEdge *)n->data;
        if (e != this)
        {
            nx = e->oppositeNode(n1);
            if (nx->mask) {nx->edges.removeNode(e); e->makeUnlinked();}
            nx->mask = 1;
        }
    }

    n = n1->edges.head();
    while (n != nullptr)
    {
        e = (GraphEdge *)n->data;
        n = n->next();
        if (e->isUnlinked()) n1->edges.removeCell((n!=nullptr)?(n->prev()):n1->edges.tail());
    }

    FOR_EACH_NODE(n1->edges, n)
            ((GraphEdge *)n->data)->oppositeNode(n1)->mask = 0;

    n1->edges.removeNode(this);

    makeUnlinked();
}


Graph::~Graph()
{
    GraphNode *gn;
    GraphEdge *ge;
    while ((gn=(GraphNode *)nodes.popHead())!=nullptr) delete gn;
    while ((ge=(GraphEdge *)edges.popHead())!=nullptr) delete ge;
}

void Graph::deleteUnlinkedElements()
{
    Node *n;
    GraphNode *gn;
    GraphEdge *ge;

    n = nodes.head();
    while (n != nullptr)
    {
        gn = (GraphNode *)n->data;
        n = n->next();
        if (gn->isIsolated())
        {
            nodes.removeCell((n!=nullptr)?(n->prev()):nodes.tail());
            delete(gn);
        }
    }

    n = edges.head();
    while (n != nullptr)
    {
        ge = (GraphEdge *)n->data;
        n = n->next();
        if (ge->isUnlinked())
        {
            edges.removeCell((n!=nullptr)?(n->prev()):edges.tail());
            delete(ge);
        }
    }
}

void Graph::unlinkEdge(GraphEdge *e)
{
    e->n1->edges.removeNode(e);
    e->n2->edges.removeNode(e);
    e->makeUnlinked();
}

void Graph::destroyEdge(GraphEdge *e)
{
    unlinkEdge(e);
    edges.removeNode(e);
    delete e;
}

GraphNode *Graph::unlinkNode(GraphNode *a)
{
    GraphEdge *e;
    while ((e = (GraphEdge *)a->edges.popHead())!=nullptr) unlinkEdge(e);
    return a;
}

void GraphEdge::invert()
{
    GraphNode *tmp = n1;
    n1 = n2;
    n2 = tmp;
}

bool Graph::isConnected()
{
    if (nodes.numberElements() < 2) return true;

    unsigned char *nmask = new unsigned char[nodes.numberElements()];
    Node *n;
    GraphNode *p, *q;
    int i;

    for (i=0, n=nodes.head(); n!=nullptr; n=n->next(), i++)
    {
        p = (GraphNode *)n->data;
        nmask[i]=p->mask;
        p->mask=0;
    }

    p = (GraphNode *)nodes.head()->data;
    List todo(p); p->mask = 1;
    while ((p = (GraphNode *)todo.popHead())!=nullptr)
    {
        for (n=p->edges.head(); n!=nullptr; n=n->next())
        {
            q = ((GraphEdge *)n->data)->oppositeNode(p);
            if (q->mask==0) {todo.appendTail(q); q->mask=1;}
        }
    }

    bool is_connected = true;
    for (i=0, n=nodes.head(); n!=nullptr; n=n->next(), i++)
    {
        p = (GraphNode *)n->data;
        if (p->mask==0) is_connected=false;
        p->mask = nmask[i];
    }

    return is_connected;
}

} 
