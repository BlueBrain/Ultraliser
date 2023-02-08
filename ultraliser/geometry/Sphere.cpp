#include "Sphere.h"
#include <iostream>


float SQUARE(float NUM) {return (NUM * NUM);}



namespace Ultraliser
{
Sphere::Sphere(Vector3f center, float radius)
{
    _center = center;
    _radius = radius;
}

bool Sphere::intersectsLine(Vector3f p1, Vector3f p2)
{
    const float a = SQUARE(p2.x() - p1.x()) + SQUARE(p2.y() - p1.y()) + SQUARE(p2.z() - p1.z());

    const float b = 2 * ((p2.x() - p1.x()) * (p1.x() - _center.x()) +
                         (p2.y() - p1.y()) * (p1.y() - _center.y()) +
                         (p2.z() - p1.z()) * (p1.z() - _center.z()));

    const float c = (SQUARE(_center.x()) + SQUARE(_center.y()) + SQUARE(_center.z()) +
                     SQUARE(p1.x()) + SQUARE(p1.y()) + SQUARE(p1.z()) -
                     2 * ((_center.x() * p1.x()) + (_center.y() * p1.y()) + (_center.z() * p1.z())) -
                     SQUARE(_radius));

    const float val = (b * b) - (4 * a * c);
    if (val < 0.f)
        return false;

    return true;
}


bool Sphere::intersectsTriangle(Vector3f p1, Vector3f p2, Vector3f p3)
{
    if (intersectsLine(p1, p2))
        return true;

    if (intersectsLine(p2, p3))
        return true;

    if (intersectsLine(p3, p1))
        return true;

    return false;
}

}
