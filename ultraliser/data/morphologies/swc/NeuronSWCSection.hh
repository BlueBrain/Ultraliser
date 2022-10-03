#pragma

#include <data/morphologies/swc/NeuronSWCSample.hh>

namespace Ultraliser
{

/**
 * @brief The NeuronSWCSection class
 */
struct NeuronSWCSection
{
    /**
     * @brief samples
     * A list of SWC samples.
     */
    NeuronSWCSamples samples;

    /**
     * @brief index
     * The index of the section.
     */
    size_t index;

    /**
     * @brief parentIndex
     * The index of the parent SWC section
     */
    size_t parentIndex;

    /**
     * @brief childrenIndices
     * A list of the indices of the children.
     */
    std::vector< size_t > childrenIndices;

    ~NeuronSWCSection()
    {
        samples.clear();
        samples.shrink_to_fit();
        childrenIndices.clear();
        childrenIndices.shrink_to_fit();
    }
};

/**
 * @brief NeuronSWCSections
 * A list of SWC sections.
 */
typedef std::vector< NeuronSWCSection* > NeuronSWCSections;

}
