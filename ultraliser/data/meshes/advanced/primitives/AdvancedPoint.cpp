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

#include "AdvancedPoint.h"
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <data/meshes/advanced/Defines.hh>


namespace Ultraliser
{

/**
 * @brief INFINITE_POINT
 */
const AdvancedPoint INFINITE_POINT(DBL_MAX, DBL_MAX, DBL_MAX);

double orient2D(const double& px, const double& py,
                const double& qx, const double& qy,
                const double& rx, const double& ry)
{
    double pqr[6];

    pqr[0] = double(px);  pqr[1] = double(py);
    pqr[2] = double(qx);  pqr[3] = double(qy);
    pqr[4] = double(rx);  pqr[5] = double(ry);

    return orient2d(pqr, pqr + 2, pqr + 4);
}

double orient3D(const AdvancedPoint *t,
                const AdvancedPoint *a,
                const AdvancedPoint *b,
                const AdvancedPoint *c)
{
    double p1[3], p2[3], p3[3], p4[3];

    p1[0] = double(t->x); p1[1] = double(t->y); p1[2] = double(t->z);
    p2[0] = double(a->x); p2[1] = double(a->y); p2[2] = double(a->z);
    p3[0] = double(b->x); p3[1] = double(b->y); p3[2] = double(b->z);
    p4[0] = double(c->x); p4[1] = double(c->y); p4[2] = double(c->z);

    return orient3d(p1, p2, p3, p4);
}

double AdvancedPoint::exactOrientation(const AdvancedPoint *a,
                                       const AdvancedPoint *b,
                                       const AdvancedPoint *c) const
{
    return orient3D(this, a, b, c);
}

bool AdvancedPoint::exactMisalignment(const AdvancedPoint *A, const AdvancedPoint *B) const
{
    if (orient2D(x, y, A->x, A->y, B->x, B->y) != 0) return true;
    if (orient2D(y, z, A->y, A->z, B->y, B->z) != 0) return true;
    if (orient2D(z, x, A->z, A->x, B->z, B->x) != 0) return true;

    return false;
}

bool AdvancedPoint::exactSameSideOnPlane(const AdvancedPoint *Q,
                                         const AdvancedPoint *A,
                                         const AdvancedPoint *B) const
{
    double o1, o2;
    int s1, s2;

    o1 = orient2D(x, y, A->x, A->y, B->x, B->y);
    o2 = orient2D(Q->x, Q->y, A->x, A->y, B->x, B->y);
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));

    if (s1 != s2) return false;

    o1 = orient2D(y, z, A->y, A->z, B->y, B->z);
    o2 = orient2D(Q->y, Q->z, A->y, A->z, B->y, B->z);
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));

    if (s1 != s2) return false;

    o1 = orient2D(z, x, A->z, A->x, B->z, B->x);
    o2 = orient2D(Q->z, Q->x, A->z, A->x, B->z, B->x);
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    return true;
}

// Returns true if the coplanar point 'p' is in the inner area of 't'.
// Undetermined if p and t are not coplanar.
bool AdvancedPoint::pointInInnerTriangle(const AdvancedPoint *p,
                                         const AdvancedPoint *v1,
                                         const AdvancedPoint *v2,
                                         const AdvancedPoint *v3)
{
    // if (!p->exactSameSideOnPlane(v1, v2, v3)) return false;
    // if (!p->exactSameSideOnPlane(v2, v3, v1)) return false;
    // if (!p->exactSameSideOnPlane(v3, v1, v2)) return false;
    // return true;

    // The following code is less readable, but slightly more efficient (12 predicates instead of 18)

    double o1, o2, oo2, oo4, oo6;
    int s1, s2;

    o1 = orient2D(p->x, p->y, v2->x, v2->y, v3->x, v3->y);
    o2 = oo2 = orient2D(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y);
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->y, p->z, v2->y, v2->z, v3->y, v3->z);
    o2 = oo4 = orient2D(v1->y, v1->z, v2->y, v2->z, v3->y, v3->z);
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->z, p->x, v2->z, v2->x, v3->z, v3->x);
    o2 = oo6 = orient2D(v1->z, v1->x, v2->z, v2->x, v3->z, v3->x);
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->x, p->y, v3->x, v3->y, v1->x, v1->y);
    o2 = oo2, s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->y, p->z, v3->y, v3->z, v1->y, v1->z);
    o2 = oo4;
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->z, p->x, v3->z, v3->x, v1->z, v1->x);
    o2 = oo6;
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->x, p->y, v1->x, v1->y, v2->x, v2->y);
    o2 = oo2;
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->y, p->z, v1->y, v1->z, v2->y, v2->z);
    o2 = oo4;
    s1 = (o1 > 0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    o1 = orient2D(p->z, p->x, v1->z, v1->x, v2->z, v2->x);
    o2 = oo6;
    s1 = (o1 > 0) ? (1) : ((o1 < 0) ? (-1) : (0)); s2 = (o2 > 0) ? (1) : ((o2 < 0) ? (-1) : (0));
    if (s1 != s2) return false;

    return true;
}

bool AdvancedPoint::operator<(const AdvancedPoint& s) const
{
    if (x < s.x)
        return true;
    else if (x>s.x)
        return false;

    if (y < s.y)
        return true;
    else if (y>s.y)
        return false;

    if (z < s.z)
        return true;
    else
        return false;
}

int xyzCompare(const void *a, const void *b)
{
    double c;

    if ((c=(((AdvancedPoint *)a)->x - ((AdvancedPoint *)b)->x)) < 0)
    {
        return -1;
    }

    if (c > 0)
    {
        return 1;
    }

    if ((c=(((AdvancedPoint *)a)->y - ((AdvancedPoint *)b)->y)) < 0)
    {
        return -1;
    }

    if (c > 0)
    {
        return 1;
    }

    if ((c=(((AdvancedPoint *)a)->z - ((AdvancedPoint *)b)->z)) < 0)
    {
        return -1;
    }

    if (c > 0)
    {
        return 1;
    }

    return 0;
}

void AdvancedPoint::normalize()
{
    double l = length();

    if (l == 0)
    {
        LOG_ERROR("AdvancedPoint::normalize : Trying to normalize a nullptr vector!");
    }

    x /= l;
    y /= l;
    z /= l;
}

void AdvancedPoint::rotate(const AdvancedPoint& a, const double& ang)
{
    double l, q[4], m[3][3];
    if ((l = a.length())==0.0) return;
    l = sin(ang/2.0)/l;

    q[0] = double(a.x) * l;
    q[1] = double(a.y) * l;
    q[2] = double(a.z) * l;
    q[3] = cos(ang/2.0);

    m[0][0] = 1.0 - (q[1] * q[1] + q[2] * q[2]) * 2.0;
    m[0][1] = (q[0] * q[1] + q[2] * q[3]) * 2.0;
    m[0][2] = (q[2] * q[0] - q[1] * q[3]) * 2.0;

    m[1][0] = (q[0] * q[1] - q[2] * q[3]) * 2.0;
    m[1][1] = 1.0 - (q[2] * q[2] + q[0] * q[0]) * 2.0;
    m[1][2] = (q[1] * q[2] + q[0] * q[3]) * 2.0;

    m[2][0] = (q[2] * q[0] + q[1] * q[3]) * 2.0;
    m[2][1] = (q[1] * q[2] - q[0] * q[3]) * 2.0;
    m[2][2] = 1.0 - (q[1] * q[1] + q[0] * q[0]) * 2.0;

    q[0] = double(x); q[1] = double(y); q[2] = double(z);
    x = m[0][0] * q[0] + m[1][0] * q[1] + m[2][0] * q[2];
    y = m[0][1] * q[0] + m[1][1] * q[1] + m[2][1] * q[2];
    z = m[0][2] * q[0] + m[1][2] * q[1] + m[2][2] * q[2];
}

void AdvancedPoint::project(const AdvancedPoint *nor)
{
    AdvancedPoint pr = (*this) - ((*nor) * ((*this) * (*nor)));
    x = pr.x;
    y = pr.y;
    z = pr.z;
}

double AdvancedPoint::distanceFromLine(const AdvancedPoint *A, const AdvancedPoint *B) const
{
    AdvancedPoint BA = (*B) - (*A);
    double lba = BA.length();

    if (lba == 0.0)
    {
        LOG_ERROR("AdvancedPoint::distanceFromLine : Degenerate line passed !");
    }

    return ((((*this) - (*A))&BA).length()) / (lba);
}

double AdvancedPoint::distanceFromLine(const AdvancedPoint *A,
                                       const AdvancedPoint *B,
                                       AdvancedPoint *cc) const
{
    AdvancedPoint AB = (*A) - (*B);
    AdvancedPoint AP = (*A) - (*this);
    AdvancedPoint BP = (*B) - (*this);

    if (AP.isNull())
    {
        cc->x = A->x;
        cc->y = A->y;
        cc->z = A->z;
        return 0.0;
    }
    else if (BP.isNull())
    {
        cc->x = B->x;
        cc->y = B->y;
        cc->z = B->z;
        return 0.0;
    }

    double t = (AB * AB);
    if (t == 0.0)
    {
        LOG_ERROR("AdvancedPoint::distanceFromLine : Degenerate line passed !");
    }
    else t = (AP * AB) / (-t);
    cc->x = t * AB.x + A->x;
    cc->y = t * AB.y + A->y;
    cc->z = t * AB.z + A->z;
    return distanceFromLine(A, B);
}

AdvancedPoint AdvancedPoint::projection(const AdvancedPoint *A, const AdvancedPoint *B) const
{
    AdvancedPoint BA = (*B)-(*A);
    double l = BA * BA;

    if (l == 0.0)
    {
        LOG_ERROR("projection : Degenerate line passed !");
    }

    return ((*A)+(BA*((BA*((*this)-(*A)))/(l))));
}

double AdvancedPoint::distanceFromEdge(const AdvancedPoint *A, const AdvancedPoint *B) const
{
    AdvancedPoint AP = (*A)-(*this);
    double apl = AP.length();
    AdvancedPoint BP = (*B) - (*this);
    double bpl = BP.length();

    if (apl == 0 || bpl == 0.0)
    {
        return 0.0;
    }

    AdvancedPoint AB = (*A) - (*B);
    double abl = AP.length();
    AdvancedPoint BA = (*B)-(*A);

    if (abl * apl == 0.0 || abl * bpl == 0.0)
    {
        return apl;
    }

    if (AB.getAngle(AP) > PI2) return apl;
    else if (BA.getAngle(BP) > PI2) return bpl;

    return distanceFromLine(A,B);
}

double AdvancedPoint::distanceFromEdge(const AdvancedPoint *A,
                                       const AdvancedPoint *B,
                                       AdvancedPoint *cc) const
{
    AdvancedPoint AP = (*A)-(*this);
    double apl = AP.length();
    AdvancedPoint BP = (*B) - (*this);
    double bpl = BP.length();

    if (apl == 0)
    {
        cc->setValue(A);
        return 0.0;
    }

    if (bpl == 0)
    {
        cc->setValue(B);
        return 0.0;
    }

    AdvancedPoint AB = (*A)-(*B); double abl = AP.length();
    AdvancedPoint BA = (*B)-(*A);

    if (abl * apl == 0.0 || abl * bpl == 0.0)
    {
        cc->setValue(A);
        return apl;
    }

    if (AB.getAngle(AP) > PI2)
    {
        cc->setValue(A);
        return apl;
    }
    else if (BA.getAngle(BP) > PI2)
    {
        cc->setValue(B);
        return bpl;
    }

    double t = (AB * AB);
    if (t == 0.0)
    {
        cc->setValue(A);
        return apl;
    }
    else
    {
        t = (AP * AB) / (-t);
    }

    cc->x = t * AB.x + A->x;
    cc->y = t * AB.y + A->y;
    cc->z = t * AB.z + A->z;
    return distanceFromLine(A, B);
}

double AdvancedPoint::getAngle(const AdvancedPoint& p) const
{
    return atan2(((*this)&p).length(), double(((*this) * p)));
}

double AdvancedPoint::distanceLineLine(const AdvancedPoint *A,
                                       const AdvancedPoint *A1,
                                       const AdvancedPoint *B1) const
{
    AdvancedPoint uu1 = ((*this)-(*A)) & ((*A1)-(*B1));
    double nom = ((*A)-(*A1)) * (uu1);
    return FABS(double(nom)) / (uu1.length());
}

AdvancedPoint AdvancedPoint::linearSystem(const AdvancedPoint& a,
                                          const AdvancedPoint& b,
                                          const AdvancedPoint& c)
{
    AdvancedPoint ret;
    const double detA = DETERMINANT3X3(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
    if (detA == 0.0)
    {
        return INFINITE_POINT;
    }

    ret.x = DETERMINANT3X3(x, a.y, a.z, y, b.y, b.z, z, c.y, c.z);
    ret.y = DETERMINANT3X3(a.x, x, a.z, b.x, y, b.z, c.x, z, c.z);
    ret.z = DETERMINANT3X3(a.x, a.y, x, b.x, b.y, y, c.x, c.y, z);

    return (ret/detA);
}

int AdvancedPoint::closestPoints(const AdvancedPoint *v1,
                                 const AdvancedPoint *p1,
                                 const AdvancedPoint *p2,
                                 AdvancedPoint *ptOnThis,
                                 AdvancedPoint *ptOnLine2) const
{
    AdvancedPoint pos1 = *this;
    AdvancedPoint dir1 = (*v1) - pos1;
    AdvancedPoint pos2 = *p1;
    AdvancedPoint dir2 = (*p2) - pos2;

    const double d1l = dir1.length();
    const double d2l = dir2.length();

    if (d1l == 0.0 && d2l == 0.0)
    {
        ptOnThis->setValue(this);
        ptOnLine2->setValue(p1);
        return 1;
    }

    if (d1l * d2l == 0.0)
    {
        if (d1l <= d2l)
        {
            ptOnThis->setValue(this);
            distanceFromLine(p1, p2, ptOnLine2);
            return 1;
        }

        if (d2l <= d1l)
        {
            ptOnLine2->setValue(p1);
            p1->distanceFromLine(this, v1, ptOnThis);
            return 1;
        }
    }

    double ang = dir1.getAngle(dir2);
    if (ang == 0.0 || ang == M_PI)
    {
        return 0;
    }

    double s, t, A, B, C, D, E, F, denom;
    denom = ((dir1 * dir2)/(d1l * d2l));
    denom = denom * denom - 1;

    dir1.normalize();
    dir2.normalize();

    A = E = dir1*dir2;
    B = dir1*dir1;
    C = (dir1*pos1) - (dir1*pos2);
    D = dir2*dir2;
    F = (dir2*pos1) - (dir2*pos2);

    s = ( C * D - A * F ) / denom;
    t = ( C * E - B * F ) / denom;
    *ptOnThis  = pos1 + (dir1 * s);
    *ptOnLine2 = pos2 + (dir2 * t);

    // Uncomment the following to compute the distance between segments
    // if (s < 0 || s > ((*v1)-(*this)).length() || t < 0 || t > ((*p2)-(*p1)).length())
    // The points does not belong to the edges
    // return 0;

    return 1;
}

AdvancedPoint AdvancedPoint::lineLineIntersection(const AdvancedPoint& p,
                                                  const AdvancedPoint& q,
                                                  const AdvancedPoint& r,
                                                  const AdvancedPoint& s)
{
    const AdvancedPoint da = q - p;
    const AdvancedPoint db = s - r;
    const AdvancedPoint dc = r - p;
    const AdvancedPoint dab = (da&db);

    if (dc * dab != 0.0)
        return INFINITE_POINT;

    double k = (((dc & db) * dab) / (dab * dab));
    return p + (da * k);
}

AdvancedPoint AdvancedPoint::linePlaneIntersection(const AdvancedPoint& p,
                                                   const AdvancedPoint& q,
                                                   const AdvancedPoint& r,
                                                   const AdvancedPoint& s,
                                                   const AdvancedPoint& t)
{
    const double denominator = DETERMINANT3X3(p.x - q.x, p.y - q.y,  p.z - q.z,
                                              s.x - r.x, s.y - r.y, s.z - r.z,
                                              t.x - r.x, t.y - r.y, t.z - r.z);

    if (denominator == 0)
        return INFINITE_POINT;

    const double numerator = DETERMINANT3X3(p.x - r.x, p.y - r.y, p.z - r.z,
                                            s.x - r.x, s.y - r.y, s.z - r.z,
                                            t.x - r.x, t.y - r.y, t.z - r.z);
    const double gamma = numerator / denominator;
    return p + ((q - p) * gamma);
}

double AdvancedPoint::squaredTriangleArea3D(const AdvancedPoint& p,
                                            const AdvancedPoint& q,
                                            const AdvancedPoint& r)
{
    AdvancedPoint pr = (p - r), qr = (q - r);
    AdvancedPoint n = pr & qr;
    return static_cast<double>((n * n) / 4.f);
}

bool AdvancedPoint::pointInInnerSegment(const AdvancedPoint *p,
                                        const AdvancedPoint *v1,
                                        const AdvancedPoint *v2)
{
    // Segment and point aligned
    if (!p->exactMisalignment(v1, v2))
    {
        if (v1->x < v2->x && v1->x < p->x && p->x < v2->x)
        {
            return true;
        }
        if (v1->y < v2->y && v1->y < p->y && p->y < v2->y)
        {
            return true;
        }
        if (v1->z < v2->z && v1->z < p->z && p->z < v2->z)
        {
            return true;
        }
        if (v1->x > v2->x && v1->x > p->x && p->x > v2->x)
        {
            return true;
        }
        if (v1->y > v2->y && v1->y > p->y && p->y > v2->y)
        {
            return true;
        }
        if (v1->z > v2->z && v1->z > p->z && p->z > v2->z)
        {
            return true;
        }
    }
    return false;
}

bool AdvancedPoint::pointInSegment(const AdvancedPoint *p,
                                   const AdvancedPoint *v1,
                                   const AdvancedPoint *v2)
{
    return ((*p) == (*(v1)) || (*p) == (*(v2)) || AdvancedPoint::pointInInnerSegment(p, v1, v2));
}

bool AdvancedPoint::pointInTriangle(const AdvancedPoint *p,
                                    const AdvancedPoint *v1,
                                    const AdvancedPoint *v2,
                                    const AdvancedPoint *v3)
{
    if (AdvancedPoint::pointInSegment(p, v1, v2))
    {
        return true;
    }
    else if (AdvancedPoint::pointInSegment(p, v2, v3))
    {
        return true;
    }
    else if (AdvancedPoint::pointInSegment(p, v3, v1))
    {
        return true;
    }
    else
    {
        return AdvancedPoint::pointInInnerTriangle(p, v1, v2, v3);
    }
}

bool AdvancedPoint::segmentsIntersect(const AdvancedPoint *p1,
                                      const AdvancedPoint *p2,
                                      const AdvancedPoint *sp1,
                                      const AdvancedPoint *sp2)
{
    return (p1->exactOrientation(p2, sp1, sp2) == 0 &&
            !p1->exactSameSideOnPlane(p2, sp1, sp2) && !sp1->exactSameSideOnPlane(sp2, p1, p2));
}

bool AdvancedPoint::innerSegmentsCross(const AdvancedPoint& p1,
                                       const AdvancedPoint& p2,
                                       const AdvancedPoint& sp1,
                                       const AdvancedPoint& sp2)
{
    if (p1 == sp1 || p1 == sp2 || p2 == sp1 || p2 == sp2)
    {
        return false;
    }

    return (p1.exactOrientation(&p2, &sp1, &sp2) == 0 &&
            !p1.exactSameSideOnPlane(&p2, &sp1, &sp2) && !sp1.exactSameSideOnPlane(&sp2, &p1, &p2));
}

bool AdvancedPoint::segmentIntersectsTriangle(const AdvancedPoint *s1,
                                              const AdvancedPoint *s2,
                                              const AdvancedPoint *v1,
                                              const AdvancedPoint *v2,
                                              const AdvancedPoint *v3)
{
    double minX = MIN(s1->x, s2->x);
    if (v1->x < minX && v2->x < minX && v3->x < minX)
    {
        return false;
    }

    minX = MAX(s1->x, s2->x);
    if (v1->x > minX && v2->x > minX && v3->x > minX)
    {
        return false;
    }

    minX = MIN(s1->y, s2->y);
    if (v1->y < minX && v2->y < minX && v3->y < minX)
    {
        return false;
    }

    minX = MAX(s1->y, s2->y);
    if (v1->y > minX && v2->y > minX && v3->y > minX)
    {
        return false;
    }

    minX = MIN(s1->z, s2->z);
    if (v1->z < minX && v2->z < minX && v3->z < minX)
    {
        return false;
    }

    minX = MAX(s1->z, s2->z);
    if (v1->z > minX && v2->z > minX && v3->z > minX)
    {
        return false;
    }

    double o1 = s1->exactOrientation(v1, v2, v3);
    double o2 = s2->exactOrientation(v1, v2, v3);

    if (o1 == 0 && o2 == 0)
    {
        if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2))
        {
            return true;
        }

        if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2))
        {
            return true;
        }

        if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2))
        {
            return true;
        }

        if (AdvancedPoint::pointInInnerTriangle(s1, v1, v2, v3) &&
            AdvancedPoint::pointInInnerTriangle(s2, v1, v2, v3))
        {
            return true;
        }

        return false;
    }

    if ((o1>0 && o2>0) || (o1<0 && o2<0))
    {
        return false;
    }

    o1 = s1->exactOrientation(s2, v1, v2);
    o2 = s1->exactOrientation(s2, v2, v3);

    if ((o1>0 && o2<0) || (o1<0 && o2>0))
    {
        return false;
    }

    double o3 = s1->exactOrientation(s2, v3, v1);
    if ((o1>0 && o3<0) || (o1<0 && o3>0))
    {
        return false;
    }

    if ((o2>0 && o3<0) || (o2<0 && o3>0))
    {
        return false;
    }

    return true;
}

bool AdvancedPoint::segmentIntersectsTriangle(const AdvancedPoint *s1,
                                              const AdvancedPoint *s2,
                                              const AdvancedPoint *v1,
                                              const AdvancedPoint *v2,
                                              const AdvancedPoint *v3,
                                              const double& oo1, const double& oo2)
{
    if (oo1 == 0 && oo2 == 0)
    {
        if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2))
        {
            return true;
        }

        if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2))
        {
            return true;
        }

        if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2))
        {
            return true;
        }

        if (AdvancedPoint::pointInInnerTriangle(s1, v1, v2, v3) &&
            AdvancedPoint::pointInInnerTriangle(s2, v1, v2, v3))
        {
            return true;
        }

        return false;
    }

    if ((oo1 > 0 && oo2 > 0) || (oo1 < 0 && oo2 < 0))
    {
        return false;
    }

    const double o1 = s1->exactOrientation(s2, v1, v2);
    const double o2 = s1->exactOrientation(s2, v2, v3);

    if ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0))
        return false;

    const double o3 = s1->exactOrientation(s2, v3, v1);

    if ((o1 > 0 && o3 < 0) || (o1 < 0 && o3 > 0))
    {
        return false;
    }

    if ((o2 > 0 && o3 < 0) || (o2 < 0 && o3 > 0))
    {
        return false;
    }

    return true;
}

}
