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

#include "AdvancedVertex.h"
#include "AdvancedEdge.h"
#include "AdvancedTriangle.h"
#include <common/Headers.hh>

namespace Ultraliser
{

AdvancedVertex::AdvancedVertex() : AdvancedPoint()
{
    e0 = nullptr;
    mask = 0;
}

AdvancedVertex::AdvancedVertex(const double &a, const double &b, const double &c)
    : AdvancedPoint(a, b, c)
{
    e0 = nullptr;
    mask = 0;
}

AdvancedVertex::AdvancedVertex(const AdvancedPoint *p)
    : AdvancedPoint(p->x, p->y, p->z)
{
    e0 = nullptr;
    mask = 0;
}

AdvancedVertex::AdvancedVertex(const AdvancedPoint& p)
    : AdvancedPoint(p.x, p.y, p.z)
{
    e0 = nullptr;
    mask = 0;
}

AdvancedVertex::~AdvancedVertex()
{
    /// EMPTY
}

/**
 * @brief Vertex::VE
 * This function should only be used for MANIFOLD and ORIENTED meshes,
 * otherwise it will never return.
 * @return
 */
List *AdvancedVertex::getIncidentEdges() const
{
    AdvancedTriangle *triangle;
    AdvancedEdge *edge;
    AdvancedVertex *vertex;
    List *incidentEdges = new List();

    if (e0 == nullptr)
        return incidentEdges;

    edge = e0;
    do
    {
        incidentEdges->appendTail(edge);
        vertex = edge->oppositeVertex(this);
        triangle = edge->leftTriangle(this);

        if (triangle == nullptr)
            break;

        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    if (edge == e0 && incidentEdges->numberElements() > 1)
        return incidentEdges;

    incidentEdges->popHead();
    edge = e0;

    do
    {
        incidentEdges->appendHead(edge);
        vertex = edge->oppositeVertex(this);
        triangle = edge->rightTriangle(this);

        if (triangle == nullptr)
            break;

        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    return incidentEdges;
}

List *AdvancedVertex::getAdjacentVertices() const
{
    AdvancedTriangle *triangle;
    AdvancedEdge *edge;
    AdvancedVertex *vertex;
    List *adjacentVertices = new List();

    if (e0 == nullptr)
        return adjacentVertices;

    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        adjacentVertices->appendTail(vertex);
        triangle = edge->leftTriangle(this);

        if (triangle == nullptr)
            break;

        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    if (edge == e0 && adjacentVertices->numberElements() > 1)
        return adjacentVertices;

    adjacentVertices->popHead();
    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        adjacentVertices->appendHead(vertex);
        triangle = edge->rightTriangle(this);

        if (triangle == nullptr)
            break;

        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    return adjacentVertices;
}

List *AdvancedVertex::VT() const
{
    AdvancedTriangle *triangle;
    AdvancedEdge *edge;
    AdvancedVertex *vertex;
    List *vt = new List();

    if (e0 == nullptr)
        return vt;

    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        triangle = edge->leftTriangle(this);

        if (triangle == nullptr)
            break;

        vt->appendTail(triangle);
        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    if (edge == e0 && vt->numberElements() > 1)
        return vt;

    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        triangle = edge->rightTriangle(this);

        if (triangle == nullptr)
            break;

        vt->appendHead(triangle);
        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    return vt;
}

AdvancedEdge *AdvancedVertex::getEdge(const AdvancedVertex *v2) const
{
    List *incidentEdges = getIncidentEdges();
    Node *node;
    AdvancedEdge *edge;

    FOR_EACH_VE_EDGE(incidentEdges, edge, node)
    {
        if (edge->oppositeVertex(this) == v2)
        {
            delete(incidentEdges);
            return edge;
        }
    }

    delete(incidentEdges);
    return nullptr;
}

int AdvancedVertex::valence() const
{
    List *incidentEdges = getIncidentEdges();

    int numberElements = incidentEdges->numberElements();

    delete(incidentEdges);

    return numberElements;
}

int AdvancedVertex::isOnBoundary() const
{
    AdvancedTriangle *triangles;
    AdvancedEdge *edge;
    AdvancedVertex *vertex;

    if (e0 == nullptr)
        return 0;

    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        triangles = edge->leftTriangle(this);

        if (triangles == nullptr)
            return 1;

        edge = triangles->oppositeEdge(vertex);
    } while (edge != e0);

    return 0;
}

AdvancedEdge *AdvancedVertex::nextBoundaryEdge() const
{
    AdvancedTriangle *triangle;
    AdvancedEdge *edge;
    AdvancedVertex *vertex;

    if (e0 == nullptr)
        return nullptr;

    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        triangle = edge->leftTriangle(this);

        if (triangle == nullptr)
            return edge;

        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    return nullptr;
}

AdvancedVertex *AdvancedVertex::nextOnBoundary() const
{
    AdvancedEdge *nextEdge = nextBoundaryEdge();

    if (nextEdge != nullptr)
        return nextEdge->oppositeVertex(this);

    return nullptr;
}

AdvancedEdge *AdvancedVertex::prevBoundaryEdge() const
{
    AdvancedTriangle *triangle;
    AdvancedEdge *edge;
    AdvancedVertex *vertex;

    if (e0 == nullptr)
        return nullptr;

    edge = e0;
    do
    {
        vertex = edge->oppositeVertex(this);
        triangle = edge->rightTriangle(this);

        if (triangle == nullptr)
            return edge;

        edge = triangle->oppositeEdge(vertex);
    } while (edge != e0);

    return nullptr;
}

AdvancedVertex *AdvancedVertex::prevOnBoundary() const
{
    AdvancedEdge *boundaryEdge = prevBoundaryEdge();

    if (boundaryEdge != nullptr)
        return boundaryEdge->oppositeVertex(this);

    return nullptr;
}

bool AdvancedVertex::isFlat() const
{
    List *incidentEdges = getIncidentEdges();
    Node *nodes;
    AdvancedEdge *edges;

    FOR_EACH_VE_EDGE(incidentEdges, edges, nodes)
    {
        if (edges->getConvexity() != 0)
        {
            delete incidentEdges;
            return false;
        }
    }

    delete incidentEdges;
    return true;
}

bool AdvancedVertex::isDoubleFlat(AdvancedEdge **e1, AdvancedEdge **e2) const
{
    List *incidentEdges = getIncidentEdges();
    Node *node;
    AdvancedEdge *edge;


    // Initially
    *e1 = *e2 = nullptr;

    int numberEdges = 0;
    FOR_EACH_VE_EDGE(incidentEdges, edge, node)
    {
        if (edge->getConvexity() != 0)
        {
            if (++numberEdges > 2)
            {
                delete incidentEdges;
                return false;
            }
            else if (numberEdges == 1)
            {
                *e1 = edge;
            }
            else
            {
                *e2 = edge;
            }
        }
    }
    delete incidentEdges;

    // This means that vertex is flat
    if (numberEdges == 0)
        return true;

    // This should not be possible, but just in case
    if (numberEdges == 1)
        return false;

    return (!((*e1)->oppositeVertex(this)->exactMisalignment(this, (*e2)->oppositeVertex(this))));
}

bool AdvancedVertex::removeIfRedundant(bool checkNeighborhood)
{
    AdvancedEdge *e, *e1 = nullptr, *e2 = nullptr;
    if (!isDoubleFlat(&e1, &e2)) return false;

    Node *n;
    AdvancedVertex *vo1, *vo2;
    List *ve;

    if (checkNeighborhood)
    {
        ve = VT();
        AdvancedTriangle *t;
        FOR_EACH_VT_TRIANGLE(ve, t, n)
        {
            if (t->isExactlyDegenerate())
            {
                delete ve;
                return false;
            }
        }

        delete ve;
        ve = getIncidentEdges();

        FOR_EACH_VE_EDGE(ve, e, n)
        {
            if (e->overlaps())
            {
                delete ve;
                return false;
            }
        }

    }
    else
        ve = getIncidentEdges();

    // isDoubleFlat()
    if (e1 != nullptr)
        ve->removeNode(e1);

    // isDoubleFlat()
    if (e2 != nullptr) ve->removeNode(e2);
    if (e2 != nullptr && *e1->oppositeVertex(this) == *e2->oppositeVertex(this))
    {
        delete ve;
        return false;
    }

    while (1)
    {
        FOR_EACH_VE_EDGE(ve, e, n)
        {
            vo1 = e->t1->oppositeVertex(e);
            vo2 = e->t2->oppositeVertex(e);
            if (!e->v1->exactSameSideOnPlane(e->v2, vo1, vo2) &&
                vo1->exactMisalignment(e->v1, vo2) &&
                vo1->exactMisalignment(e->v2, vo2))
            {
                if (!e->swap())
                {
                    delete ve;
                    return false;
                }
                break;
            }
        }

        if (n != nullptr)
            ve->removeCell(n);
        else
            break;
    }

    // isDoubleFlat()
    if (e1 != nullptr)
        e = e1;
    else
        if (ve->numberElements() == 3)
            e = (AdvancedEdge *)ve->head()->data;
    else
    {
        FOR_EACH_VE_EDGE(ve, e, n)
        {
            if (!e->t1->oppositeVertex(e)->exactMisalignment(
                        this, e->t2->oppositeVertex(e))) break;
        }

        if (n == nullptr)
        {
            delete ve;
            return false;
        }
        else
            e = (AdvancedEdge *)((n == ve->head()) ? (ve->tail()) : (n->prev()))->data;
    }
    delete ve;

    if (e->v1 == this)
        e->invert();
    if (e->collapseOnV1() == nullptr)
        return false;
    else return true;
}

AdvancedPoint AdvancedVertex::getNormal() const
{
    List *vertices = VT();
    Node *node;
    AdvancedTriangle *triangle;
    AdvancedPoint tnor, triangleNormal;

    double angle;

    FOR_EACH_VT_TRIANGLE(vertices, triangle, node)
    {
        angle = triangle->getAngle(this);
        triangleNormal = triangle->getNormal();

        if (!triangleNormal.isNull())
            tnor = tnor+(triangleNormal * angle);
    }
    delete(vertices);

    if (tnor.isNull())
        return AdvancedPoint(0,0,0);

    tnor.normalize();
    return tnor;
}

double AdvancedVertex::getBoundaryAngle() const
{
    AdvancedEdge *e1 = prevBoundaryEdge();
    AdvancedEdge *e2 = nextBoundaryEdge();

    if (e1 == nullptr || e2 == nullptr)
        return -1.0;

    AdvancedVertex *v1 = e1->oppositeVertex(this);
    AdvancedVertex *v2 = e2->oppositeVertex(this);
    const double angle = getAngle(v1, v2);

    return angle;
}

double AdvancedVertex::getAngleForTriangulation() const
{
    AdvancedEdge *e1 = prevBoundaryEdge();
    AdvancedEdge *e2 = nextBoundaryEdge();

    if (e1 == nullptr || e2 == nullptr)
        return DBL_MAX;

    AdvancedTriangle *t1 = e1->getBoundaryTriangle();
    AdvancedTriangle *t2 = e2->getBoundaryTriangle();

    AdvancedVertex *v1 = e1->oppositeVertex(this);
    AdvancedVertex *v2 = e2->oppositeVertex(this);

    if ((*v2)==(*v1))
        return -2;

    if (distance(v1)*distance(v2) == 0.0)
        return -1;

    double angle = getAngle(v1, v2);
    if (angle == M_PI)
        return 3 * M_PI;
    if (angle == 0)
        return 0;

    AdvancedEdge e3(v1,v2);
    AdvancedTriangle t(e1,e2,&e3);
    double dihAngle1 = t.getDAngle(t1);
    double dihAngle2 = t.getDAngle(t2);

    if (dihAngle1==M_PI && dihAngle2==M_PI)
        return (DBL_MAX * 0.5f);

    if (dihAngle1==M_PI || dihAngle2==M_PI)
        return (DBL_MAX / 4.0f);

    return dihAngle1+dihAngle2+angle;
}

double AdvancedVertex::getAngleOnAveragePlane(AdvancedPoint *normal) const
{
    AdvancedEdge *e1 = prevBoundaryEdge();
    AdvancedEdge *e2 = nextBoundaryEdge();

    if (e1 == nullptr || e2 == nullptr)
        return DBL_MAX;

    AdvancedVertex *v1 = e1->oppositeVertex(this);
    AdvancedVertex *v2 = e2->oppositeVertex(this);

    AdvancedPoint p, p1, p2;
    p1.setValue(v1);
    p2.setValue(v2);
    p.setValue(this);
    p.project(normal);
    p1.project(normal);
    p2.project(normal);

    if (p.distance(p1)*p.distance(p2) == 0.0)
    {
        LOG_WARNING("getAngleOnAveragePlane: coincident projections!");
        return 0.0;
    }

    double angle = p.getAngle(&p1, &p2);
    if (normal->side3D(&p1, &p, &p2) < 0)
        angle = -(angle - (2 * M_PI));

    return angle;
}

double AdvancedVertex::totalDihedralAngle() const
{
    List *incidentTriangles = getIncidentEdges();
    double mc = 0;
    AdvancedEdge *edge;
    Node *node;

    FOR_EACH_VE_EDGE(incidentTriangles, edge, node)
    {
        if (edge->isOnBoundary())
        {
            delete(incidentTriangles);
            return DBL_MAX;
        }
        else
        {
            mc -= (edge->dihedralAngle() - M_PI);
        }
    }
    mc /= incidentTriangles->numberElements();

    delete(incidentTriangles);

    return mc;
}

double AdvancedVertex::totalAngle() const
{
    List *incidentTriangles = getIncidentEdges();
    double totalAngleF = 0.0;

    AdvancedEdge *iEdge;
    Node *node;

    FOR_EACH_VE_EDGE(incidentTriangles, iEdge, node)
    {
        if (iEdge->isOnBoundary())
        {
            delete(incidentTriangles);
            return -1.0;
        }
        else
        {
            totalAngleF += iEdge->leftTriangle(this)->getAngle(this);
        }
    }
    delete(incidentTriangles);

    return totalAngleF;
}

double AdvancedVertex::voronoiArea() const
{
    List *incidentTriangles = VT();
    Node *node;
    AdvancedTriangle *triangle;

    double totalvoronoiArea = 0.0;
    FOR_EACH_VT_TRIANGLE(incidentTriangles, triangle, node)
    {
        totalvoronoiArea += triangle->area();
    }
    delete(incidentTriangles);

    return totalvoronoiArea / 3.0f;
}

int AdvancedVertex::zip(const bool checkGeometry)
{
    Node *node;
    AdvancedEdge *iEdge;

    List *incidentEdges = getIncidentEdges();

    AdvancedEdge *be1 = (AdvancedEdge *) incidentEdges->head()->data;
    AdvancedEdge *be2 = (AdvancedEdge *) incidentEdges->tail()->data;
    delete(incidentEdges);

    if (!be1->isOnBoundary() || !be2->isOnBoundary())
        return 0;

    AdvancedVertex *ov1 = be1->oppositeVertex(this);
    AdvancedVertex *ov2 = be2->oppositeVertex(this);

    if (checkGeometry && ((*ov1)!=(*ov2)))
        return 0;

    if (ov1 != ov2)
    {
        incidentEdges = ov2->getIncidentEdges();

        FOR_EACH_VE_EDGE(incidentEdges, iEdge, node)
        {
            iEdge->replaceVertex(ov2, ov1);
        }

        delete(incidentEdges);
        ov2->e0 = nullptr;
    }

    AdvancedTriangle *triangle = (be2->t1!=nullptr) ? (be2->t1) : (be2->t2);
    triangle->replaceEdge(be2, be1);
    be1->replaceTriangle(nullptr, triangle);

    be2->v1 = be2->v2 = nullptr;

    e0 = ov1->e0 = be1;

    return 1 + ov1->zip(checkGeometry);
}

AdvancedEdge *AdvancedVertex::inverseCollapse(AdvancedVertex *v2,
                                              AdvancedEdge *e,
                                              AdvancedEdge *e1,
                                              AdvancedEdge *e2,
                                              AdvancedEdge *e3,
                                              AdvancedEdge *e4,
                                              AdvancedTriangle *t1,
                                              AdvancedTriangle *t2)
{
    AdvancedTriangle *rightTriangle = e2->rightTriangle(this);
    AdvancedTriangle *leftTriangle = e3->leftTriangle(this);

    Node *node;
    AdvancedEdge *iEdge;
    List *incidentEdges = getIncidentEdges();
    FOR_EACH_VE_EDGE(incidentEdges, iEdge, node)
    {
        if (iEdge == e3)
        {
            break;
        }
    }

    FOR_EACH_NODE_CIRCULAR((*incidentEdges), node, node)
    {
        iEdge = ((AdvancedEdge *) node->data);

        if (iEdge == e2)
            break;
        else
            iEdge->replaceVertex(this, v2);
    }
    delete(incidentEdges);

    e->v1 = this;
    e->v2 = v2;

    e1->v1 = v2;
    e1->v2 = e2->oppositeVertex(this);
    e4->v1 = v2;
    e4->v2 = e3->oppositeVertex(this);

    t1->edge1 = e;
    t1->edge2 = e1;
    t1->edge3 = e2;
    t2->edge1 = e;
    t2->edge2 = e3;
    t2->edge3 = e4;

    e->t1 = t1;
    e->t2 = t2;
    e2->replaceTriangle(rightTriangle, t1);
    e3->replaceTriangle(leftTriangle, t2);

    if (rightTriangle)
        rightTriangle->replaceEdge(e2, e1);
    if (leftTriangle)
        leftTriangle->replaceEdge(e3, e4);

    e1->t1 = t1;
    e1->t2 = rightTriangle;
    e4->t1 = leftTriangle;
    e4->t2 = t2;
    v2->e0 = e0 = e;

    return e;
}

}
