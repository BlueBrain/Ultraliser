#pragma once

#include <data/meshes/simple/Mesh.h>

namespace Ultraliser
{

class IcoSphere : public Mesh
{
public:
    IcoSphere(const size_t& subdivisions);
    ~IcoSphere() { };

    size_t getNumberSubdivisions() const { return _subdivisions; };

public:

    size_t _subdivisions;
};

}
