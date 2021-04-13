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

#include "AdvancedTriangle.h"
#include <data/meshes/advanced/Defines.hh>


namespace Ultraliser
{

extern "C" double orient2d(double *, double *, double *);

AdvancedTriangle::AdvancedTriangle(){
    mask = 0;
    info = nullptr;
}

AdvancedTriangle::AdvancedTriangle(AdvancedEdge *a, AdvancedEdge *b, AdvancedEdge *c)
{
    edge1 = a;
    edge2 = b;
    edge3 = c;
    mask = 0;
    info = nullptr;
}

AdvancedPoint AdvancedTriangle::getNormal() const
{
    AdvancedVertex *va = v1();
    AdvancedVertex *vb = v2();
    AdvancedVertex *vc = v3();
    AdvancedPoint vd = (((*va)-(*vb)) & ((*vb)-(*vc)));

    const double l = vd.length();

    if (l == 0)
        return AdvancedPoint(0.f, 0.f, 0.f);

    return vd / l;
}

AdvancedPoint AdvancedTriangle::getVector() const
{
    AdvancedVertex *va = v1();
    AdvancedVertex *vb = v2();
    AdvancedVertex *vc = v3();

    return (((*va) - (*vb)) & ((*vb) - (*vc)));
}

bool AdvancedTriangle::checkAdjNor(const AdvancedTriangle *t) const
{
    AdvancedEdge *e = commonEdge(t);
    if (e == nullptr)
        return true;

    AdvancedEdge *ea = nextEdge(e);
    AdvancedEdge *eb = t->nextEdge(e);

    if (ea->commonVertex(eb) == ea->commonVertex(e))
        return false;

    return true;
}

double AdvancedTriangle::area() const
{
    const double a = edge1->length();
    const double b = edge2->length();
    const double c = edge3->length();

    if (a == 0.0 || b == 0.0 || c == 0.0)
        return 0.0;

    double p = (a + b + c) * 0.5f;

    p = p * (p - a) * (p - b) * (p - c);

    if (p < 0)
        return 0.0;

    return sqrt(p);
}

double AdvancedTriangle::perimeter() const
{
    return edge1->length() + edge2->length() + edge3->length();
}

AdvancedPoint AdvancedTriangle::getCenter() const
{
    AdvancedPoint va = *v1(), vb = *v2(), vc = *v3();
    return (va + vb + vc) / 3.0f;
}

AdvancedPoint AdvancedTriangle::getCircleCenter() const
{
    AdvancedPoint va = *v1(), vb = *v2(), vc = *v3();
    AdvancedPoint q1 = vb-va;
    AdvancedPoint q2 = vc-va;
    AdvancedPoint n = q2 & q1;
    AdvancedPoint m1 = edge2->getMidPoint();
    AdvancedPoint m2 = edge1->getMidPoint();

    return AdvancedPoint(n * va,q1 * m1,q2 * m2).linearSystem(n, q1, q2);
}

bool AdvancedTriangle::inSphere(const AdvancedPoint *p) const
{
    AdvancedPoint center = getCircleCenter();
    double radius = center.squaredDistance(edge1->v1);

    return (p->squaredDistance(&center) < radius);
}

double AdvancedTriangle::getAngle(const AdvancedVertex *v) const
{
    AdvancedVertex *va = v1();
    AdvancedVertex *vb = v2();
    AdvancedVertex *vc = v3();

    if (v == va)
        return v->getAngle(vb, vc);

    if (v == vb)
        return v->getAngle(va, vc);

    if (v == vc)
        return v->getAngle(vb, va);

    return -1.0;
}

double AdvancedTriangle::getDAngle(const AdvancedTriangle *t) const
{
    AdvancedPoint thisNormal = getVector();
    AdvancedPoint otherNormal = t->getVector();

    if (thisNormal.isNull() || otherNormal.isNull())
        return -1.0;

    return thisNormal.getAngle(otherNormal);
}

double AdvancedTriangle::distanceFromPoint(const AdvancedPoint *p) const
{
    return sqrt(double(squaredDistanceFromPoint(p)));
}

double AdvancedTriangle::squaredDistanceFromPoint(const AdvancedPoint *p) const 
{
    AdvancedPoint CA = edge1->toVector() & edge2->toVector();
    double CA2 = CA * CA;

    if (CA2 == 0)
        return -1.0;

    const double d = ((CA * (*p))-(CA * (*(edge1->v1))));

    return (d * d) / CA2;
}

double AdvancedTriangle::pointTriangleDistance(const AdvancedPoint *p, AdvancedPoint *cp) const
{
    return sqrt(double(pointTriangleSquaredDistance(p)));
}

double AdvancedTriangle::pointTriangleSquaredDistance(const AdvancedPoint *p,
                                                      AdvancedEdge **closestEdge,
                                                      AdvancedVertex **closestVertex) const
{
    AdvancedVertex *va = v1(), *vb = v2(), *vc = v3();
    AdvancedPoint n(((*va) - (*vb)) & ((*vb) - (*vc)));

    if (n.x == 0 && n.y == 0 && n.z == 0)
        return -1.0;

    double d1 = ((((*va)-(*vb)) & ((*vb)-(*p))) * n);
    double d2 = ((((*vb)-(*vc)) & ((*vc)-(*p))) * n);
    double d3 = ((((*vc)-(*va)) & ((*va)-(*p))) * n);

    // Closest point in inner triangle
    if (d1 > 0 && d2 > 0 && d3 > 0)
    {
        if (closestEdge != nullptr)
            *closestEdge = nullptr;

        if (closestVertex != nullptr)
            *closestVertex = nullptr;

        return squaredDistanceFromPoint(p);
    }

    if (d2 < 0)
    {
        va = vb;
        vb = vc;

        if (closestEdge != nullptr)
            *closestEdge = edge3;
    }
    else if (d3 < 0)
    {
        vb = va;
        va = vc;

        if (closestEdge != nullptr)
            *closestEdge = edge1;
    }
    else if (closestEdge != nullptr)
        *closestEdge = edge2;

    AdvancedPoint i(p->projection(va, vb));
    AdvancedPoint p1(i - (*va));
    AdvancedPoint p2(i - (*vb));

    // Closest point on interior of one edge
    if (p1 * p2 < 0)
    {
        return i.squaredDistance(p);
    }

    d1 = p1.squaredLength();
    d2 = p2.squaredLength();

    if (d1 < d2)
    {
        if (closestVertex != nullptr)
            *closestVertex = va;
        return p->squaredDistance(va);
    }
    else
    {
        if (closestVertex != nullptr)
            *closestVertex = vb;
        return p->squaredDistance(vb);
    }
}

AdvancedPoint AdvancedTriangle::project(const AdvancedPoint *p) const
{
    AdvancedPoint n = getVector();

    if (n.isNull())
        return INFINITE_POINT;

    return AdvancedPoint::linePlaneIntersection(*p, (*p) + n, *(v1()), *(v2()), *(v3()));
}

bool AdvancedTriangle::isExactlyDegenerate() const
{
    return (!v1()->exactMisalignment((v2()), (v3())));
}

AdvancedEdge *AdvancedTriangle::getLongestEdge() const
{
    const double l1 = edge1->squaredLength();
    const double l2 = edge2->squaredLength();
    const double l3 = edge3->squaredLength();

    if (l1 >= l2 && l1 >= l3)
        return edge1;

    if (l2 >= l1 && l2 >= l3)
        return edge2;

    return edge3;
}

AdvancedEdge *AdvancedTriangle::getCapEdge() const
{
    AdvancedEdge *edge;

    edge = edge1;
    if (AdvancedPoint::pointInInnerSegment(oppositeVertex(edge), edge->v1, edge->v2))
        return edge;

    edge = edge2;
    if (AdvancedPoint::pointInInnerSegment(oppositeVertex(edge), edge->v1, edge->v2))
        return edge;

    edge = edge3;
    if (AdvancedPoint::pointInInnerSegment(oppositeVertex(edge), edge->v1, edge->v2))
        return edge;

    return nullptr;
}

bool AdvancedTriangle::overlaps() const
{
    return (edge1->overlaps() || edge2->overlaps() || edge3->overlaps());
}

AdvancedVertex *AdvancedTriangle::commonVertex(const AdvancedTriangle *t2) const
{
    if (hasVertex(t2->v1()))
        return t2->v1();

    if (hasVertex(t2->v2()))
        return t2->v2();

    if (hasVertex(t2->v3()))
        return t2->v3();

    return nullptr;
}

void AdvancedTriangle::printTriangleToFile(FILE *file) const
{
    v1()->printPointToFile(file);
    v2()->printPointToFile(file);
    v3()->printPointToFile(file);
}

bool AdvancedTriangle::intersects(const AdvancedTriangle *t2, bool justProper) const
{
    AdvancedVertex *v11, *v12, *v13, *v21, *v22, *v23;

    if (justProper)
    {
        // This works for non-degenerate triangles. Not sure it will work for degeneracies too.
        v11 = v1();
        v12 = v2();
        v13 = v3();

        v21 = t2->v1();
        v22 = t2->v2();
        v23 = t2->v3();

        AdvancedVertex *equation1 = ((*v11) == (*v21)) ? (v21) :
                                   (((*v11) == (*v22)) ? (v22) :
                                   (((*v11) == (*v23)) ? (v23) : (nullptr)));
        AdvancedVertex *equation2 = ((*v12) == (*v21)) ? (v21) :
                                   (((*v12) == (*v22)) ? (v22) :
                                   (((*v12) == (*v23)) ? (v23) : (nullptr)));
        AdvancedVertex *equation3 = ((*v13) == (*v21)) ? (v21) :
                                   (((*v13) == (*v22)) ? (v22) :
                                   (((*v13) == (*v23)) ? (v23) : (nullptr)));

        // Triangles coincide
        if (equation1 && equation2 && equation3)
        {
            return false;
        }

        AdvancedEdge *ce1 = nullptr, *ce2 = nullptr;
        if (equation1 && equation2)
        {
            ce1 = edge2;
            ce2 = (t2->edge1->hasVertices(equation1, equation2)) ? (t2->edge1) :
                 ((t2->edge2->hasVertices(equation1, equation2)) ? (t2->edge2) : (t2->edge3));
        }

        if (equation2 && equation3)
        {
            ce1 = edge3;
            ce2 = (t2->edge1->hasVertices(equation3, equation2)) ? (t2->edge1) :
                 ((t2->edge2->hasVertices(equation3, equation2)) ? (t2->edge2) : (t2->edge3));
        }

        if (equation3 && equation1)
        {
            ce1 = edge1;
            ce2 = (t2->edge1->hasVertices(equation3, equation1)) ? (t2->edge1) :
                 ((t2->edge2->hasVertices(equation3, equation1)) ? (t2->edge2) : (t2->edge3));
        }

        if (ce1)
        {
            AdvancedVertex *ov = t2->oppositeVertex(ce2);

            return (ov->exactOrientation(v11, v12, v13) == 0 &&
                    ov->exactSameSideOnPlane(oppositeVertex(ce1), ce1->v1, ce1->v2));
        }

        AdvancedVertex *cv1 = nullptr, *cv2 = nullptr;
        if (equation1)
        {
            cv1 = v11;
            cv2 = equation1;
        }

        if (equation2)
        {
            cv1 = v12;
            cv2 = equation2;
        }

        if (equation3)
        {
            cv1 = v13;
            cv2 = equation3;
        }

        // If they share a vertex, intersection occurs if the opposite edge intersect the triangle
        if (cv1)
        {
            AdvancedEdge *ee1 = oppositeEdge(cv1);
            AdvancedEdge *ee2 = t2->oppositeEdge(cv2);

            return (AdvancedPoint::segmentIntersectsTriangle(ee1->v1, ee1->v2, v21, v22, v23) ||
                    AdvancedPoint::segmentIntersectsTriangle(ee2->v1, ee2->v2, v11, v12, v13));
        }
    }
    else
    {
        AdvancedEdge *ce = commonEdge(t2);
        // If they share an edge, intersection occurs only if t1 and t2 overlap
        if (ce)
        {
            AdvancedVertex *ov = t2->oppositeVertex(ce);
            return (ov->exactOrientation(v1(), v2(), v3()) == 0 &&
                    ov->exactSameSideOnPlane(oppositeVertex(ce), ce->v1, ce->v2));
        }

        AdvancedVertex *cv = commonVertex(t2);
        v11 = v1();
        v12 = v2();
        v13 = v3();

        v21 = t2->v1();
        v22 = t2->v2();
        v23 = t2->v3();

        // If they share a vertex, intersection occurs if the opposite edge intersect the triangle
        if (cv)
        {
            AdvancedEdge *ee1 = oppositeEdge(cv);
            AdvancedEdge *ee2 = t2->oppositeEdge(cv);

            return (AdvancedPoint::segmentIntersectsTriangle(ee1->v1, ee1->v2, v21, v22, v23) ||
                    AdvancedPoint::segmentIntersectsTriangle(ee2->v1, ee2->v2, v11, v12, v13));
        }
    }

    // Fast reject by bounding box
    double minValue = MIN(v11->x, MIN(v13->x, v12->x));
    if (v21->x < minValue && v22->x < minValue && v23->x < minValue)
    {
        return false;
    }

    minValue = MAX(v11->x, MAX(v13->x, v12->x));
    if (v21->x > minValue && v22->x > minValue && v23->x > minValue)
    {
        return false;
    }

    minValue = MIN(v11->y, MIN(v13->y, v12->y));
    if (v21->y < minValue && v22->y < minValue && v23->y < minValue)
    {
        return false;
    }

    minValue = MAX(v11->y, MAX(v13->y, v12->y));
    if (v21->y > minValue && v22->y > minValue && v23->y > minValue)
    {
        return false;
    }

    minValue = MIN(v11->z, MIN(v13->z, v12->z));
    if (v21->z < minValue && v22->z < minValue && v23->z < minValue)
    {
        return false;
    }

    minValue = MAX(v11->z, MAX(v13->z, v12->z));
    if (v21->z > minValue && v22->z > minValue && v23->z > minValue)
    {
        return false;
    }

    // Calculate relative orientations
    const double o11 = v11->exactOrientation(v21, v22, v23);
    const double o12 = v12->exactOrientation(v21, v22, v23);
    const double o13 = v13->exactOrientation(v21, v22, v23);

    // t1 above / below t2
    if ((o11 > 0 && o12 > 0 && o13 > 0) || (o11 < 0 && o12 < 0 && o13 < 0))
        return false;

    const double o21 = v21->exactOrientation(v11, v12, v13);
    const double o22 = v22->exactOrientation(v11, v12, v13);
    const double o23 = v23->exactOrientation(v11, v12, v13);

    // t2 above / below t1
    if ((o21 > 0 && o22 > 0 && o23 > 0) || (o21 < 0 && o22 < 0 && o23 < 0))
        return false;

    // t1 and t2 are coplanar
    if (o11 == 0 && o12 == 0 && o13 == 0)
    {
        if (AdvancedPoint::innerSegmentsCross(v11, v12, v21, v22))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v11, v12, v22, v23))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v11, v12, v23, v21))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v12, v13, v21, v22))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v12, v13, v22, v23))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v12, v13, v23, v21))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v13, v11, v21, v22))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v13, v11, v22, v23))
            return true;
        if (AdvancedPoint::innerSegmentsCross(v13, v11, v23, v21))
            return true;

        return (AdvancedPoint::pointInTriangle(v11, v21, v22, v23) ||
                AdvancedPoint::pointInTriangle(v12, v21, v22, v23) ||
                AdvancedPoint::pointInTriangle(v13, v21, v22, v23) ||
                AdvancedPoint::pointInTriangle(v21, v11, v12, v13) ||
                AdvancedPoint::pointInTriangle(v22, v11, v12, v13) ||
                AdvancedPoint::pointInTriangle(v23, v11, v12, v13));
    }
    else
        return (AdvancedPoint::segmentIntersectsTriangle(v11, v12, v21, v22, v23, o11, o12) ||
                AdvancedPoint::segmentIntersectsTriangle(v12, v13, v21, v22, v23, o12, o13) ||
                AdvancedPoint::segmentIntersectsTriangle(v13, v11, v21, v22, v23, o13, o11) ||
                AdvancedPoint::segmentIntersectsTriangle(v21, v22, v11, v12, v13, o21, o22) ||
                AdvancedPoint::segmentIntersectsTriangle(v22, v23, v11, v12, v13, o22, o23) ||
                AdvancedPoint::segmentIntersectsTriangle(v23, v21, v11, v12, v13, o23, o21));
}

}
