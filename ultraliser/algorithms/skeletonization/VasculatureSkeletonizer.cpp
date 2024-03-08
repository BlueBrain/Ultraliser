#include "VasculatureSkeletonizer.h"

namespace Ultraliser
{

VasculatureSkeletonizer::VasculatureSkeletonizer(Volume* volume, const bool &useAcceleration)
    : Skeletonizer(volume, useAcceleration)

{
    /// EMPTY CONSTRUCTOR
}

void VasculatureSkeletonizer::skeletonizeVolumeToCenterLines()
{
    // Start the timer
    TIMER_SET;

    _applyVolumeThinning();
    constructGraph();
    segmentComponents();

    LOG_STATUS_IMPORTANT("Skeletonization Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void VasculatureSkeletonizer::segmentComponents()
{
    // Start the timer
    TIMER_SET;

    LOG_STATUS("Segmenting Graph Components");

    // Currently, vasculature datasets have only branches with fluctuating diameters. Therefore,
    // branches are only the available components to be segmented from the graph. This is indeed
    // unlike neurons that can have somata, branches and spines.
    _buildBranchesFromNodes(_nodes);

     LOG_STATS(GET_TIME_SECONDS);
}

void VasculatureSkeletonizer::exportSkeletonVMV(const std::string& filePrefix)
{
    // Construct the file path
    std::string filePath = filePrefix + VMV_EXTENSION;
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

void VasculatureSkeletonizer::exportSkeletonVMV5(const std::string& filePrefix)
{
    // Construct the file path
    std::string filePath = filePrefix + VMV5_EXTENSION;
    LOG_STATUS("Exporting VMV5 Morphology : [ %s ]", filePath.c_str());

    // Open the file
    H5::H5File file(filePath, H5F_ACC_TRUNC);

    // Add the vertices dataset
    const int RANK = 2;
    hsize_t  verticesDimensions[2];
    verticesDimensions[0] = _nodes.size();
    verticesDimensions[1] = 4; // X Y Z R
    H5::DataSpace dataspace(RANK, verticesDimensions);

    // Define the datatype for the data in the file.
    H5::FloatType datatype(H5::PredType::NATIVE_FLOAT);

    // Create a new dataset container to store the dataset in.
    H5::DataSet dataset = file.createDataSet("vertices", datatype, dataspace );

    struct NodeVert {
        float x;
        float y;
        float z;
        float r;
    };

    std::vector< NodeVert > verts;
    verts.resize(_nodes.size());
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        verts[i].x = _nodes[i]->point.x();
        verts[i].y = _nodes[i]->point.y();
        verts[i].z = _nodes[i]->point.z();
        verts[i].r = _nodes[i]->radius;
    }

    // Write the dataset using the default memory space
    dataset.write(verts.data(), H5::PredType::NATIVE_FLOAT);

    // TODO: Strands are missing

    // Close the file
    file.close();
}

VasculatureSkeletonizer::~VasculatureSkeletonizer()
{
    /// EMPTY
}
}
