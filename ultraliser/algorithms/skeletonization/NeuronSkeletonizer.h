#pragma once

#include <algorithms/skeletonization/Skeletonizer.h>

namespace Ultraliser
{

/**
 * @brief The NeuronSkeletonizer class
 */
class NeuronSkeletonizer : public Skeletonizer
{
public:
    NeuronSkeletonizer(const Mesh* mesh);
};
}
