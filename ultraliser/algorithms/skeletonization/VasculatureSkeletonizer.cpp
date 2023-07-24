#include "VasculatureSkeletonizer.h"

namespace Ultraliser
{

VasculatureSkeletonizer::VasculatureSkeletonizer(const Mesh *mesh, Volume* volume)
    : Skeletonizer(mesh, volume)

{

}

void VasculatureSkeletonizer::segmentComponents()
{

    // Build the branches from the nodes
    _buildBranchesFromNodes(_nodes);

    std::cout << "Branches: " << _branches.size() << "\n";
    std::cout << "Nodes (Samples): " << _nodes.size() << "\n";

    std::fstream stream;
    stream.open("/abdellah2/scratch/thinning/output/projections/branches.txt", std::ios::out);
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        stream << "start\n";

        for (auto& node: _branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        stream << "end\n";
    }
    stream.close();

    std::cout << "Branchs: " << _branches.size() << "\n";

    return;
}
void VasculatureSkeletonizer::exportSkeletonVMV(const std::string& prefix,
                                                const std::string& fileName)
{
    std::string filePath = prefix + "/" + fileName;
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
    fprintf(filePointer, "$PARAM_END\n");

    LOOP_STARTS("Writing Nodes");
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
    fprintf(filePointer, "$VERT_LIST_END\n");
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOOP_STARTS("Writing Strands");
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
}

VasculatureSkeletonizer::~VasculatureSkeletonizer()
{

}
}
