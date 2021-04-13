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

#include "AdvancedEdge.h"
#include "AdvancedTriangle.h"
#include <data/meshes/advanced/Defines.hh>

namespace Ultraliser
{


int edgeCompare(const void *a, const void *b)
{
    double lengthA = ((AdvancedEdge *) a)->squaredLength();
    double lengthB = ((AdvancedEdge *) b)->squaredLength();

    if (lengthA < lengthB)
        return -1;

    if (lengthA > lengthB)
        return 1;

    return 0;
}

int lexEdgeCompare(const void *a, const void *b)
{
    AdvancedVertex *va1 = ((AdvancedEdge *)a)->v1;
    AdvancedVertex *va2 = ((AdvancedEdge *)a)->v2;
    AdvancedVertex *vb1 = ((AdvancedEdge *)b)->v1;
    AdvancedVertex *vb2 = ((AdvancedEdge *)b)->v2;

    if (xyzCompare(va1, va2) > 0)
        p_swap((void **)&va1, (void **)&va2);

    if (xyzCompare(vb1, vb2) > 0)
        p_swap((void **)&vb1, (void **)&vb2);

    int ca = xyzCompare(va1, vb1);

    if (ca == 0)
        return xyzCompare(va2, vb2);

    return ca;
}

int vtxEdgeCompare(const void *a, const void *b)
{
    AdvancedVertex *va1 = ((AdvancedEdge *)a)->v1;
    AdvancedVertex *va2 = ((AdvancedEdge *)a)->v2;
    AdvancedVertex *vb1 = ((AdvancedEdge *)b)->v1;
    AdvancedVertex *vb2 = ((AdvancedEdge *)b)->v2;
    AdvancedVertex *tmp;

    if (va2 < va1)
    {
        tmp = va1;
        va1 = va2;
        va2 = tmp;
    }

    if (vb2 < vb1)
    {
        tmp = vb1;
        vb1 = vb2;
        vb2 = tmp;
    }

    if (va1 < vb1)
        return -1;
    if (va1 > vb1)
        return 1;
    if (va2 < vb2)
        return -1;
    if (va2 > vb2)
        return 1;
    return 0;
}

AdvancedEdge::AdvancedEdge()
{
    mask = 0;
    info = nullptr;
}

AdvancedEdge::AdvancedEdge(AdvancedVertex *va, AdvancedVertex *vb)
{
    v1 = va;
    v2 = vb;
    t1 = t2 = nullptr;
    info = nullptr;
    mask = 0;
}

AdvancedEdge::~AdvancedEdge()
{
    /// EMPTY DESTRUCTOR
}

AdvancedPoint AdvancedEdge::toUnitVector() const
{
    AdvancedPoint v = toVector();
    double l = v.length();

    if (l == 0)
    { printf("Err");
        //LOG_ERROR("Edge::toUnitVector : Degenerate Edge !\n");
    }

    return v / l;
}

AdvancedPoint AdvancedEdge::getNormal() const
{
    AdvancedPoint normal, n1, n2;

    if (t1 == nullptr || t2 == nullptr)
        return AdvancedPoint(0, 0, 0);

    n1 = t1->getNormal();
    n2 = t2->getNormal();
    normal = n1 + n2;

    if (normal.length() != 0.0)
        normal.normalize();

    return normal;
}

bool AdvancedEdge::swap(const bool fast)
{
    if (!fast && (t1 == nullptr || t2 == nullptr ||
                  t2->oppositeVertex(this)->getEdge(
                      t1->oppositeVertex(this)) != nullptr))
        return 0;

    AdvancedEdge *e1 = t1->nextEdge(this);
    AdvancedEdge *e3 = t2->nextEdge(this);

    v1->e0 = e3;
    v2->e0 = e1;

    v1 = t2->oppositeVertex(this);
    v2 = t1->oppositeVertex(this);

    t1->replaceEdge(e1, e3);
    t2->replaceEdge(e3, e1);

    t1->invert();
    t2->invert();

    e1->replaceTriangle(t1, t2);
    e3->replaceTriangle(t2, t1);

    return 1;
}

AdvancedVertex *AdvancedEdge::collapseOnV1()
{
    AdvancedEdge *e;
    Node *n;
    List *ve;
    AdvancedVertex *tv;

    AdvancedEdge *e1 = (t1 != nullptr) ? (t1->nextEdge(this)) : (nullptr);
    AdvancedEdge *e2 = (t1 != nullptr) ? (t1->prevEdge(this)) : (nullptr);
    AdvancedEdge *e3 = (t2 != nullptr) ? (t2->nextEdge(this)) : (nullptr);
    AdvancedEdge *e4 = (t2 != nullptr) ? (t2->prevEdge(this)) : (nullptr);

    AdvancedVertex *v3 = (e1 != nullptr) ? (e1->oppositeVertex(v2)) : (nullptr);
    AdvancedVertex *v4 = (e4 != nullptr) ? (e4->oppositeVertex(v2)) : (nullptr);

    AdvancedTriangle *ta1 = (e1 != nullptr) ? (e1->oppositeTriangle(t1)) : (nullptr);
    AdvancedTriangle *ta2 = (e2 != nullptr) ? (e2->oppositeTriangle(t1)) : (nullptr);
    AdvancedTriangle *ta3 = (e3 != nullptr) ? (e3->oppositeTriangle(t2)) : (nullptr);
    AdvancedTriangle *ta4 = (e4 != nullptr) ? (e4->oppositeTriangle(t2)) : (nullptr);

    if (v1->isOnBoundary() && v2->isOnBoundary())
        if (!(((ta1 || ta2) && !ta3 && !ta4) ||
              ((ta3 || ta4) && !ta1 && !ta2))) return nullptr;

    if (ta1 != nullptr && ta2 != nullptr &&
        ta1->oppositeVertex(e1) == ta2->oppositeVertex(e2))
        return nullptr;
    if (ta3 != nullptr && ta4 != nullptr &&
        ta3->oppositeVertex(e3) == ta4->oppositeVertex(e4))
        return nullptr;

    if (ta1 == nullptr && ta2 == nullptr)
        v1->e0 = e3;
    else
        v1->e0 = e2;

    if (v3 != nullptr)
        v3->e0 = e2;

    if (v4 != nullptr)
        v4->e0 = e3;

    ve = v2->getIncidentEdges();
    FOR_EACH_VE_EDGE(ve, e, n)
    {
        tv = e->oppositeVertex(v2);
        if (tv != v3 && tv != v4 &&
            tv->getEdge(v1) != nullptr)
        {
            delete(ve);
            return nullptr;
        }
    }

    FOR_EACH_VE_EDGE(ve, e, n)
    {
        if (e != this)
            e->replaceVertex(v2, v1);
    }
    delete(ve);

    if (e2 != nullptr)
        e2->replaceTriangle(t1, ta1);
    if (e3 != nullptr)
        e3->replaceTriangle(t2, ta4);

    if (ta1 != nullptr)
        ta1->replaceEdge(e1, e2);
    if (ta4 != nullptr)
        ta4->replaceEdge(e4, e3);

    // v2 must be removed
    v2->e0 = nullptr;

    // e4 must be removed
    if (e4 != nullptr) e4->v1 = e4->v2 = nullptr;

    // e1 must be removed
    if (e1 != nullptr) e1->v1 = e1->v2 = nullptr;

    // t1 must be removed
    if (t1 != nullptr) t1->edge1 = t1->edge2 = t1->edge3 = nullptr;

    // t2 must be removed
    if (t2 != nullptr) t2->edge1 = t2->edge2 = t2->edge3 = nullptr;

    if (e2 != nullptr && e2->t1 == nullptr && e2->t2 == nullptr)
    {
        v3->e0 = nullptr;
        e2->v1 = e2->v2 = nullptr;
    }

    if (e3 != nullptr && e3->t1 == nullptr && e3->t2 == nullptr)
    {
        v4->e0 = nullptr;
        e3->v1 = e3->v2 = nullptr;
    }

    // This is the remaining vertex to be returned
    v4 = v1;

    // This edge must be removed
    v2 = v1 = nullptr;

    return v4;
}

AdvancedVertex *AdvancedEdge::collapseOnV2()
{
    invert();
    return collapseOnV1();
}

bool AdvancedEdge::collapse(const AdvancedPoint& p)
{
    AdvancedVertex *r = collapseOnV1();
    if (r==nullptr)
        return false;

    // Average the collapse
    else r->setValue(&p);

    return true;
}

bool AdvancedEdge::collapse()
{
    return collapse(((*v1) + (*v2)) / 2);
}

bool AdvancedEdge::merge(AdvancedEdge *e)
{
    if (t1 && t2)
        return 0;

    if (e->t1 && e->t2)
        return 0;

    AdvancedTriangle *ot = (e->t1==nullptr) ? (e->t2) : (e->t1);

    if (ot == getBoundaryTriangle())
        return 0;

    if ((t1 && e->t1) || (t2 && e->t2))
        e->invert();

    AdvancedVertex *ov1 = e->v1, *ov2 = e->v2;
    List *ve1=nullptr, *ve2=nullptr;
    Node *n;
    AdvancedEdge *f, *f2;

    if (ov1 != v1)
    {
        ve1 = ov1->getIncidentEdges();
        FOR_EACH_VE_EDGE(ve1, f, n)
        {
            f2 = f->oppositeVertex(ov1)->getEdge(v1);
            if (f2 != nullptr && (!f2->isOnBoundary() || !f->isOnBoundary()))
            {
                delete(ve1);
                return 0;
            }
        }
    }
    if (ov2 != v2)
    {
        ve2 = ov2->getIncidentEdges();
        FOR_EACH_VE_EDGE(ve2, f, n)
        {
            f2 = f->oppositeVertex(ov2)->getEdge(v2);
            if (f2 != nullptr && (!f2->isOnBoundary() || !f->isOnBoundary()))
            {
                delete(ve1);
                delete(ve2);
                return 0;
            }
        }
    }

    if (ov1 != v1)
    {
        FOR_EACH_VE_EDGE(ve1, f, n) f->replaceVertex(ov1, v1);
        delete(ve1);
        ov1->e0 = nullptr;
    }
    if (ov2 != v2)
    {
        FOR_EACH_VE_EDGE(ve2, f, n) f->replaceVertex(ov2, v2);
        delete(ve2);
        ov2->e0 = nullptr;
    }
    ot->replaceEdge(e, this);
    ((t1==nullptr) ? (t1) : (t2)) = ot;
    v1->e0 = v2->e0 = this;
    e->v1 = e->v2 = nullptr;

    return 1;
}

double AdvancedEdge::curvature() const
{
    if (!t1 || !t2) return -1.0;
    return t1->getDAngle(t2);
}

double AdvancedEdge::dihedralAngle() const
{
    if (!t1 || !t2)
        return -1.0;

    AdvancedPoint nor1 = t1->getNormal();
    AdvancedPoint nor2 = t2->getNormal();

    if (nor1.isNull() || nor2.isNull())
        return -1.0;
    double c = nor1.getAngle(&nor2);

    AdvancedVertex *ov = t2->oppositeVertex(this);
    if (((*ov)*nor1) - ((*v1)*nor1) < 0)
        return -(c - M_PI);

    return c + M_PI;
}

double AdvancedEdge::delaunayMinAngle() const
{
    if (t1==nullptr || t2==nullptr)
        return 2 * M_PI;
    if (squaredLength()==0)
        return 0;
    if (t1->nextEdge(this)->squaredLength() == 0)
        return 0;
    if (t1->prevEdge(this)->squaredLength() == 0)
        return 0;

    double a1 = t1->getAngle(v1);
    double a2 = t1->getAngle(v2);
    double a3 = t1->getAngle(t1->oppositeVertex(this));

    if (t2->nextEdge(this)->length()==0)
        return 0;
    if (t2->prevEdge(this)->length()==0)
        return 0;

    double a4 = t2->getAngle(v1);
    double a5 = t2->getAngle(v2);
    double a6 = t2->getAngle(t2->oppositeVertex(this));

    if (a1 + a4 >= M_PI || a2 + a5 >= M_PI)
        return 3 * M_PI;
    return MIN(a1, (MIN(a2, (MIN(a3, (MIN(a4, (MIN(a5, a6)))))))));
}

bool AdvancedEdge::stitch()
{
    // This function seems to be insufficient to stitch in every case !
    if (!isOnBoundary()) return 0;

    AdvancedTriangle *t, *t0 = (t1 != nullptr) ? (t1) : (t2);
    AdvancedVertex *v0;
    AdvancedEdge *e1;

    for (v0 = v1; v0 != nullptr; v0 = ((v0 == v1) ? (v2) : (nullptr)))
    {
        e1 = this;
        t = t0;
        while (t != nullptr)
        {
            e1 = t->nextEdge(e1);
            if (!e1->hasVertex(v0))
                e1 = t->nextEdge(e1);
            t = e1->oppositeTriangle(t);
        }
        if (e1->oppositeVertex(v0) == oppositeVertex(v0))
        {
            t = (e1->t1 != nullptr) ? (e1->t1) : (e1->t2);
            t->replaceEdge(e1, this);
            v1->e0 = v2->e0 = this;
            e1->v1 = e1->v2 = nullptr;
            replaceTriangle(nullptr, t);
            return 1;
        }
    }

    return 0;
}

bool AdvancedEdge::overlaps() const
{
    if (t1 == nullptr || t2 == nullptr)
        return false;

    AdvancedVertex *ov = t2->oppositeVertex(this);
    if (ov->exactOrientation(t1->v1(), t1->v2(), t1->v3()) == 0 &&
        ov->exactSameSideOnPlane(t1->oppositeVertex(this), v1, v2))
        return true;
    else
        return false;
}

bool AdvancedEdge::intersects(const AdvancedTriangle *t) const
{
    if (t->hasEdge(this))
        return false;

    AdvancedVertex *cv = (t->hasVertex(v1)) ? (v1) :
                ((t->hasVertex(v2)) ? (v2) : (nullptr));

    // If they share a vertex, intersection occurs if t's opposite edge
    // intersect this edge
    if (cv)
    {
        AdvancedEdge *oe = t->oppositeEdge(cv);
        if (AdvancedPoint::pointInTriangle(oppositeVertex(cv), cv, oe->v1, oe->v2))
            return true;
        else
            return (AdvancedPoint::segmentsIntersect(oe->v1, oe->v2, v1, v2));
    }
    else return AdvancedPoint::segmentIntersectsTriangle(v1, v2, t->v1(), t->v2(), t->v3());
}

double AdvancedEdge::getConvexity() const
{
    if (t1 == nullptr || t2 == nullptr)
        return DBL_MAX;
    else
        return (t1->oppositeVertex(this)->exactOrientation(t2->v3(), t2->v2(), t2->v1()));
}

}
