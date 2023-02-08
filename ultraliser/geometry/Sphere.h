#pragma once

#include <math/Vector3f.h>

namespace Ultraliser
{
class Sphere
{
public:
    Sphere(Vector3f center, float radius);

    Vector3f getCenter()  const {return _center; }

    float getRadius() const { return _radius;}

    bool intersectsLine(Vector3f p1, Vector3f p2);

    bool intersectsTriangle(Vector3f p1, Vector3f p2, Vector3f p3);


private:

    Vector3f _center;

    float _radius;
};

}
