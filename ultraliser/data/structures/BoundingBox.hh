#pragma once

#include <math/Vector3f.h>
#include <common/Headers.hh>

namespace Ultraliser
{

struct BoundingBox
{
public:

    BoundingBox(const Vector3f& pMin, const Vector3f& pMax)
    {
        this->pMin = pMin;
        this->pMax = pMax;
        this->bounds = pMax - pMin;
        this->center = pMin + bounds * 0.5;
    }


    std::vector< Vector3f > getPoint() const
    {
        std::vector< Vector3f > points;

        points.push_back(pMin);
        return points;
    }

public:

    Vector3f pMin;
    Vector3f pMax;
    Vector3f bounds;
    Vector3f center;

    std::vector< Vector3f > points;

};

}
