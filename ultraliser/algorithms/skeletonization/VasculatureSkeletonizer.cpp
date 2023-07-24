#include "VasculatureSkeletonizer.h"

namespace Ultraliser
{

VasculatureSkeletonizer::VasculatureSkeletonizer(Volume* volume, const Mesh *mesh)
    : Skeletonizer(volume, mesh)

{
    /// EMPTY CONSTRUCTOR
}

void VasculatureSkeletonizer::skeletonizeVolume()
{
    applyVolumeThinning();
    constructGraph();
    segmentComponents();
}


void VasculatureSkeletonizer::skeletonizeVolumeBlockByBlock(const size_t& blockSize,
                                                 const size_t& numberOverlappingVoxels,
                                                 const size_t& numberZeroVoxels)
{
    thinVolumeBlockByBlock(blockSize, numberOverlappingVoxels, numberZeroVoxels);
    constructGraph();
    segmentComponents();
}

void VasculatureSkeletonizer::segmentComponents()
{
    // Currently, vasculature datasets have only branches with fluctuating diameters. Therefore,
    // branches are only the available components to be segmented from the graph. This is indeed
    // unlike neurons that can have somata, branches and spines.
    _buildBranchesFromNodes(_nodes);
}

void VasculatureSkeletonizer::exportSkeletonVMV(const std::string& prefix,
                                                const std::string& fileName)
{
    // Construct the file path
    std::string filePath = prefix + "/" + fileName + VMV_EXTENSION;
    LOG_STATUS("Exporting VMV Morphology : [ %s ]", filePath.c_str());

    // Start the timer
    TIMER_SET;

    // Open file for writing
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "w")) == nullptr)
    {
        LOG_WARNING("Cannot write [ %s ]! ", filePath.c_str());
        return;
    }

    // Header
    fprintf(filePointer, "$PARAM_BEGIN\n");
    fprintf(filePointer, "$NUM_VERTS %lu\n", _nodes.size());
    fprintf(filePointer, "$NUM_STRANDS %lu\n", _branches.size());
    fprintf(filePointer, "$PARAM_END\n\n");

    // Nodes, or vertices
    LOOP_STARTS("Writing Vertices (or Nodes)");
    fprintf(filePointer, "$VERT_LIST_BEGIN\n");
    size_t progress = 0;
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        LOOP_PROGRESS(progress, _nodes.size());
        ++progress;

        const auto& node = _nodes[i];
        const auto& index = node->index;
        const auto& p = node->point;
        const auto& r = node->radius;
        fprintf(filePointer, "%lu %f %f %f %f\n", index, p.x(), p.y(), p.z(), r);
    }
    fprintf(filePointer, "$VERT_LIST_END\n\n");
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Strands, or segments
    LOOP_STARTS("Writing Strands (or Segments)");
    fprintf(filePointer, "$STRANDS_LIST_BEGIN\n");
    progress = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        LOOP_PROGRESS(progress, _branches.size());
        ++progress;

        const auto& branch = _branches[i];
        const auto& index = branch->index;
        fprintf(filePointer, "%lu ", index);
        for (size_t j = 0; j < branch->nodes.size(); ++j)
        {
            fprintf(filePointer, "%lu ", branch->nodes[j]->index);
        }
        fprintf(filePointer, "\n");
    }
    fprintf(filePointer, "$STRANDS_LIST_END\n");
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    fclose(filePointer);

    LOG_STATUS_IMPORTANT("Writing VMV File Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

VasculatureSkeletonizer::~VasculatureSkeletonizer()
{

}
}
