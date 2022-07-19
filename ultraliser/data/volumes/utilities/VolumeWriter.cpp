#include "VolumeWriter.h"
#include <utilities/Timer.h>
#include <data/volumes/grids/Grids.h>

namespace Ultraliser
{

void writeNRRD(const std::string &prefix, const BitVolumeGrid* grid)
{
    // Starts the timer
    TIMER_SET;

    // Open the file in write mode
    std::string fileName = prefix + std::string(NRRD_EXTENSION);
    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    FILE* fptr= fopen(fileName.c_str(), "w");

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: Volume\n");
    fprintf(fptr, "type: uint16\n");
    fprintf(fptr, "dimension: 3\n");
    // fprintf(fptr, "space dimension: 3\n");
    fprintf(fptr, "sizes: %" PRId64 " %" PRId64 " %" PRId64 "\n",
            grid->getWidth(), grid->getHeight(), grid->getDepth());
    fprintf(fptr, "spacings: 1 1 1\n");
    fprintf(fptr, "endian: little\n");
    fprintf(fptr, "encoding: raw\n");

    LOOP_STARTS("Writing Voxels (1 Byte per voxel)");
    for (size_t i = 0; i < grid->getNumberVoxels(); ++i)
    {
        LOOP_PROGRESS(i, grid->getNumberVoxels());
        uint16_t value = grid->isFilled(i) ? FILLED_UI8_VOXEL_VALUE : EMPTY_VOXEL_VALUE;
        // fputc(value, fptr);
        fwrite(&value, 2, 1, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

template< class T >
void writeNRRD(const std::string &prefix, const UnsignedVolumeGrid<T>* grid)
{
    // Starts the timer
    TIMER_SET;

    // File name
    std::string fileName = prefix + std::string(NRRD_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    size_t numberBytesPerVoxel = 0;
    size_t voxelSizeInBytes = 0;

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    if (typeid (T) == typeid (uint8_t))
    {
        fprintf(fptr, "type: unsigned char\n");
        numberBytesPerVoxel = 1;
        voxelSizeInBytes = sizeof(uint8_t);
    }
    else if (typeid (T) == typeid (uint16_t))
    {
        fprintf(fptr, "type: unsigned short\n");
        numberBytesPerVoxel = 2;
        voxelSizeInBytes = sizeof(uint16_t);
    }
    else if (typeid (T) == typeid (uint32_t))
    {
        fprintf(fptr, "type: unsigned int\n");
        numberBytesPerVoxel = 4;
        voxelSizeInBytes = sizeof(uint32_t);
    }
    else if (typeid (T) == typeid (uint64_t))
    {
        fprintf(fptr, "type: unsigned long\n");
        numberBytesPerVoxel = 8;
        voxelSizeInBytes = sizeof(uint64_t);
    }
    else
    {
        LOG_ERROR("Undefined volume type!");
    }

    fprintf(fptr, "dimension: 3\n");
    fprintf(fptr, "sizes: %zu %zu %zu\n", grid->getWidth(), grid->getHeight(), grid->getDepth());
    fprintf(fptr, "spacings: 1 1 1\n");
    fprintf(fptr, "encoding: raw\n");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels");
    fwrite(grid->getGridData(), voxelSizeInBytes, grid->getNumberVoxels(), fptr);
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

template< class T >
void writeNRRD(const std::string &prefix, const FloatVolumeGrid<T>* grid)
{
    // Starts the timer
    TIMER_SET;

    // File name
    std::string fileName = prefix + std::string(NRRD_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    size_t numberBytesPerVoxel = 0;
    size_t voxelSizeInBytes = 0;

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    if (typeid (T) == typeid (float))
    {
        fprintf(fptr, "type: float\n");
        numberBytesPerVoxel = 4;
        voxelSizeInBytes = sizeof(float);
    }
    else if (typeid (T) == typeid (double))
    {
        fprintf(fptr, "type: double\n");
        numberBytesPerVoxel = 8;
        voxelSizeInBytes = sizeof(double);
    }
    else
    {
        LOG_ERROR("Undefined volume type!");
    }

    fprintf(fptr, "dimension: 3\n");
    fprintf(fptr, "sizes: %zu %zu %zu\n", grid->getWidth(), grid->getHeight(), grid->getDepth());
    fprintf(fptr, "spacings: 1 1 1\n");
    fprintf(fptr, "encoding: raw\n");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels");
    fwrite(grid->getGridData(), voxelSizeInBytes, grid->getNumberVoxels(), fptr);
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void writeVOL(const std::string &prefix, const BitVolumeGrid* grid,
              const bool &oneBitPerVoxel)
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "HEADER_BEGIN\n");
    if (oneBitPerVoxel) { fprintf(fptr, "format:%s\n", FORMAT_BIT.c_str()); }
    else { fprintf(fptr, "format:%s\n", FORMAT_8UI.c_str()); }
    fprintf(fptr, "size:%zux%zux%zu\n", grid->getWidth(), grid->getHeight(), grid->getDepth());
    fprintf(fptr, "HEADER_END\n");

    if (oneBitPerVoxel)
    {
        LOOP_STARTS("Writing Voxels (1 Bit per voxel)");
        for (size_t i = 0; i < grid->getNumberBytes(); ++i)
        {
            LOOP_PROGRESS_FRACTION(i, grid->getNumberBytes());
            fputc(grid->getByte(i), fptr);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        LOOP_STARTS("Writing Voxels (1 Byte per voxel)");
        for (size_t i = 0; i < grid->getNumberVoxels(); ++i)
        {
            LOOP_PROGRESS(i, grid->getNumberVoxels());
            uint8_t value = grid->isFilled(i) ? FILLED_UI8_VOXEL_VALUE : EMPTY_VOXEL_VALUE;
            fputc(value, fptr);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Closing the file
    fclose(fptr);
}

template< class T >
void writeVOL(const std::string &prefix,  const UnsignedVolumeGrid<T>* grid)
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "HEADER_BEGIN\n");

    size_t voxelSizeInBytes = 0;

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    if (typeid (T) == typeid (uint8_t))
    {
        fprintf(fptr, "format:%s\n", FORMAT_8UI.c_str());
        voxelSizeInBytes = sizeof(uint8_t);
    }
    else if (typeid (T) == typeid (uint16_t))
    {
        fprintf(fptr, "format:%s\n", FORMAT_16UI.c_str());
        voxelSizeInBytes = sizeof(uint16_t);
    }
    else if (typeid (T) == typeid (uint32_t))
    {
        fprintf(fptr, "format:%s\n", FORMAT_32UI.c_str());
        voxelSizeInBytes = sizeof(uint32_t);
    }
    else if (typeid (T) == typeid (uint64_t))
    {
        fprintf(fptr, "format:%s\n", FORMAT_64UI.c_str());
        voxelSizeInBytes = sizeof(uint64_t);
    }
    else
    {
        LOG_ERROR("Undefined volume type!");
    }
    fprintf(fptr, "size:%zux%zux%zu\n", grid->getWidth(), grid->getHeight(), grid->getDepth());
    fprintf(fptr, "HEADER_END\n");

    fwrite(grid->getGridData(), voxelSizeInBytes, grid->getNumberVoxels(), fptr);

    // Closing the file
    fclose(fptr);
}

template< class T >
void writeVOL(const std::string &prefix,  const FloatVolumeGrid<T>* grid)
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "HEADER_BEGIN\n");

    size_t numberBytesPerVoxel = 0;
    size_t voxelSizeInBytes = 0;

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    if (typeid (T) == typeid (float))
    {
        fprintf(fptr, "format:%s\n", FORMAT_F32.c_str());
        numberBytesPerVoxel = 4;
        voxelSizeInBytes = sizeof(float);
    }
    else if (typeid (T) == typeid (double))
    {
        fprintf(fptr, "format:%s\n", FORMAT_F64.c_str());
        numberBytesPerVoxel = 8;
        voxelSizeInBytes = sizeof(double);
    }
    else
    {
        LOG_ERROR("Undefined volume type!");
    }
    fprintf(fptr, "size:%zux%zux%zu\n", grid->getWidth(), grid->getHeight(), grid->getDepth());
    fprintf(fptr, "HEADER_END\n");

    fwrite(grid->getGridData(), voxelSizeInBytes, grid->getNumberVoxels(), fptr);

    // Closing the file
    fclose(fptr);
}

void writeHeaderFile(const std::string &prefix, const VOLUME_TYPE& type,
                     const size_t& width, const size_t& height, const size_t& depth)
{
    // Header path
    std::string fileName = prefix + std::string(HEADER_EXTENSION);

    // Open the file
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);

    LOG_STATUS("Writing Header [ %s ]", fileName.c_str());

    if (typeid (type) == typeid (uint8_t))
    {
        header << FORMAT_8UI << std::endl;
    }
    else if (typeid (type) == typeid (uint16_t))
    {
        header << FORMAT_16UI << std::endl;
    }
    else if (typeid (type) == typeid (uint32_t))
    {
        header << FORMAT_32UI << std::endl;
    }
    else if (typeid (type) == typeid (uint64_t))
    {
        header << FORMAT_64UI << std::endl;
    }

    // Write the dimensions
    header << width << " " << height << " " << depth << std::endl;

    // Close the file
    header.close();
}

void writeRAW(const std::string &prefix, const BitVolumeGrid* grid)
{
    // Starts the timer
    TIMER_SET;

    // Write the header file, always use an 8UIfor the BitGrid
    writeHeaderFile(prefix, VOLUME_TYPE::UI8,
                    grid->getWidth(), grid->getHeight(), grid->getDepth());

    std::string fileName = prefix + std::string(RAW_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Byte per voxel)");
    for (size_t i = 0; i < grid->getNumberVoxels(); ++i)
    {
        LOOP_PROGRESS(i, grid->getNumberVoxels());
        uint8_t value = grid->isFilled(i) ? FILLED_UI8_VOXEL_VALUE : EMPTY_VOXEL_VALUE;
        fputc(value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

template< class T >
void writeRAW(const std::string &prefix,  const UnsignedVolumeGrid<T>* grid)
{
    VOLUME_TYPE type;
    size_t voxelSizeInBytes = 0;

    if (typeid (T) == typeid (uint8_t))
    {
        type = VOLUME_TYPE::UI8;
        voxelSizeInBytes = sizeof(uint8_t);
    }
    else if (typeid (T) == typeid (uint16_t))
    {
        type = VOLUME_TYPE::UI16;
        voxelSizeInBytes = sizeof(uint16_t);
    }
    else if (typeid (T) == typeid (uint32_t))
    {
        type = VOLUME_TYPE::UI32;
        voxelSizeInBytes = sizeof(uint32_t);
    }
    else if (typeid (T) == typeid (uint64_t))
    {
        type = VOLUME_TYPE::UI64;
        voxelSizeInBytes = sizeof(uint64_t);
    }

    // Write the header file, always use an 8UIfor the BitGrid
    writeHeaderFile(prefix, type,
                    grid->getWidth(), grid->getHeight(), grid->getDepth());

    std::string fileName = prefix + std::string(RAW_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    fwrite(grid->getGridData(), voxelSizeInBytes, grid->getNumberVoxels(), fptr);

    // Closing the file
    fclose(fptr);
}

// Templeta Specialization
#include "VolumeWriter.ipp"

}
