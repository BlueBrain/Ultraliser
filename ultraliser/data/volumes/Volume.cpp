/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#include "Volume.h"
#include <cstddef>
#include <data/volumes/utilities/VolumeType.hh>
#include <common/Headers.hh>
#include <common/Common.h>
#include <geometry/Intersection.h>
#include <geometry/Utilities.h>
#include <utilities/Utilities.h>
#include <math/Functions.h>
#include <data/images/TIFFImage.h>
#include <data/volumes/grids/VolumeGrid.h>
#include <data/volumes/voxels/DMCVoxel.h>
#include <data/volumes/grids/BitVolumeGrid.h>
#include <data/volumes/grids/UnsignedVolumeGrid.h>
#include <data/volumes/grids/Grids.h>
#include <data/meshes/simple/VolumeMesh.h>
#include <data/meshes/simple/MeshOperations.h>
#include <data/volumes/utilities/VolumeReader.h>
#include <algorithms/skeletonization/Neighbors.hh>
#include <algorithms/skeletonization/Thinning6Iterations.h>
#include <data/meshes/simple/TriangleOperations.h>
#include <algorithm>

#ifdef ULTRALISER_USE_NRRD
#include <nrrdloader/NRRDLoader.h>
#endif

namespace Ultraliser
{

Volume::Volume(const Volume *inputVolume)
    : _surfaceVoxelizationTime(0.f)
    , _solidVoxelizationTime(0.f)
    , _addingVolumePassTime(0.f)
{
    // Copy the bounds
    _pMin = inputVolume->getPMin();
    _pMax = inputVolume->getPMax();
    _center = inputVolume->getCenter();
    _scale = inputVolume->getScale();

    // Copy the voxel size
    _voxelSize = inputVolume->getVoxelSize();

    // Copy the grid type
    _gridType = inputVolume->getType();

    switch (_gridType)
    {
    case VOLUME_TYPE::BIT:
    {
        const BitVolumeGrid* grid = static_cast< BitVolumeGrid* >(inputVolume->getGrid());
        _grid = new BitVolumeGrid(grid);
    } break;

    case VOLUME_TYPE::UI8:
    {
        const VolumeGridU8* grid = static_cast< VolumeGridU8* >(inputVolume->getGrid());
        _grid = new VolumeGridU8(inputVolume->getWidth(),
                                 inputVolume->getHeight(),
                                 inputVolume->getDepth(),
                                 grid->getGridData());
    } break;

    case VOLUME_TYPE::UI16:
    {
        const VolumeGridU16* grid = static_cast< VolumeGridU16* >(inputVolume->getGrid());
        _grid = new VolumeGridU16(inputVolume->getWidth(),
                                 inputVolume->getHeight(),
                                 inputVolume->getDepth(),
                                 grid->getGridData());
    } break;

    case VOLUME_TYPE::UI32:
    {
        const VolumeGridU32* grid = static_cast< VolumeGridU32* >(inputVolume->getGrid());
        _grid = new VolumeGridU32(inputVolume->getWidth(),
                                 inputVolume->getHeight(),
                                 inputVolume->getDepth(),
                                 grid->getGridData());
    } break;

    case VOLUME_TYPE::UI64:
    {
        const VolumeGridU64* grid = static_cast< VolumeGridU64* >(inputVolume->getGrid());
        _grid = new VolumeGridU64(inputVolume->getWidth(),
                                 inputVolume->getHeight(),
                                 inputVolume->getDepth(),
                                 grid->getGridData());
    } break;

    case VOLUME_TYPE::F32:
    {
        _gridType = VOLUME_TYPE::F32;
        const VolumeGridF32* grid = static_cast< VolumeGridF32* >(inputVolume->getGrid());
        _grid = new VolumeGridF32(inputVolume->getWidth(),
                                 inputVolume->getHeight(),
                                 inputVolume->getDepth(),
                                 grid->getGridData());

    } break;

    case VOLUME_TYPE::F64:
    {
        const VolumeGridF64* grid = static_cast< VolumeGridF64* >(inputVolume->getGrid());
        _grid = new VolumeGridF64(inputVolume->getWidth(),
                                 inputVolume->getHeight(),
                                 inputVolume->getDepth(),
                                 grid->getGridData());
    } break;

    default:
    { LOG_ERROR("Undefined volume format"); } break;
    }
}

Volume::Volume(const std::string &filePath)
    : _pMin(Vector3f::ZERO)
    , _pMax(Vector3f::ZERO)
    , _surfaceVoxelizationTime(0.f)
    , _solidVoxelizationTime(0.f)
    , _addingVolumePassTime(0.f)
{
    // Get the volume extension
    auto extension = std::filesystem::path(filePath).extension().string();

    // To lower for a single comparison
    String::toLower(extension);

    // Read the NRRD file and construct the grid
    if (String::subStringFound(extension, NRRD_EXTENSION))
    {
#ifdef ULTRALISER_USE_NRRD
        _createGrid(readNRRDVolumeFile(filePath));
#else
        LOG_ERROR("NRRD support is missing due to unavailable dependencies!");
#endif
    }

    // Read the HDR/IMG or HDR/BIN files
    else if (String::subStringFound(extension, HEADER_EXTENSION))
    {
        // Construct the grid directly from the header file
        _createGrid(filePath);
    }

    // Read the Ultraliser volume file
    else if (String::subStringFound(extension, ULTRALISER_VOLUME_EXTENSION) ||
             String::subStringFound(extension, ULTRALISER_BIN_VOLUME_EXTENSION))
    {
        // Read the volume data from the file and construct the volume grid
        _createGrid(readUltraliserVolumeFile(filePath));
    }

    // Unrecognized format
    else
    {
        LOG_ERROR("The input volume extenstion [%s] is NOT supported!", extension.c_str());
    }

    // Since we don't have any geometric bounds, use 1.0 for the voxel resolution
    // TODO: The voxel resolution should be computed from the given bounds
    _voxelSize = 1.f;

    // If the bounds are zero, then adjust the scale based on the grid dimensions.
    if (_pMin.isZero() and _pMax.isZero())
    {
        // pMin will be at the origin
        _pMin = Vector3f::ZERO;

        // pMax will be at the higher bound
        _pMax.x() = 1.f * _grid->getWidth();
        _pMax.y() = 1.f * _grid->getHeight();
        _pMax.z() = 1.f * _grid->getDepth();

        _scale = _pMax - _pMin;
        _center = _pMin + (0.5f * _scale);
    }
}

Volume::Volume(const Vector3f& pMin,
               const Vector3f& pMax,
               const size_t &baseResolution,
               const float &expansionRatio,
               const VOLUME_TYPE& gridType)
    : _gridType(gridType)
    , _pMin(pMin)
    , _pMax(pMax)
    , _expansionRatio(expansionRatio)
    , _baseResolution(baseResolution)
    , _surfaceVoxelizationTime(0.f)
    , _solidVoxelizationTime(0.f)
    , _addingVolumePassTime(0.f)
{
    // Create the grid
    _createGrid();
}

Volume::Volume(const int64_t width,
               const int64_t height,
               const int64_t depth,
               const Vector3f pMin,
               const Vector3f pMax,
               const VOLUME_TYPE& gridType,
               const float expansionRatio)
    : _gridType(gridType)
    , _pMin(pMin)
    , _pMax(pMax)
    , _expansionRatio(expansionRatio)
    , _surfaceVoxelizationTime(0.f)
    , _solidVoxelizationTime(0.f)
    , _addingVolumePassTime(0.f)
{
    // Since we don't have any geometric bounds, use 1.0 for the voxel resolution
    // TODO: The voxel resolution should be computed from the given bounds
    _voxelSize = 1.0;

    // Allocating the grid
    _allocateGrid(width, height, depth);
}

#ifdef ULTRALISER_USE_NRRD
void Volume::_createGrid(const NRRDVolumeData* volumeData)
{
    // Update the cente, scale and pMin and pMax accordingly
    _center = volumeData->center;
    _scale = volumeData->scale;
    _pMin = _center - 0.5 * _scale;
    _pMax = _center + 0.5 * _scale;

    switch (volumeData->type)
    {
    case VOLUME_TYPE::BIT:
    {
        LOG_ERROR("Unimplemented Volume::_createGrid(const NRRDVolumeData* volumeData)");
    } break;

    case VOLUME_TYPE::UI8:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI8;
        _grid = new VolumeGridU8(volumeData->width, volumeData->height, volumeData->depth,
                                 volumeData->data->asBytes().data());
    } break;

    case VOLUME_TYPE::UI16:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI16;
        _grid = new VolumeGridU16(volumeData->width, volumeData->height, volumeData->depth,
                                  volumeData->data->asUnsingedShorts().data());
    } break;

    case VOLUME_TYPE::UI32:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI32;
        _grid = new VolumeGridU32(volumeData->width, volumeData->height, volumeData->depth,
                                  volumeData->data->asUnsignedIntegers().data());
    } break;

    case VOLUME_TYPE::UI64:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI64;
        _grid = new VolumeGridU64(volumeData->width, volumeData->height, volumeData->depth,
                                  volumeData->data->asUnsignedLongs().data());
    } break;

    case VOLUME_TYPE::F32:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::F32;
        _grid = new VolumeGridF32(volumeData->width, volumeData->height, volumeData->depth,
                                  volumeData->data->asFloats().data());

    } break;

    case VOLUME_TYPE::F64:
    {
        _gridType = VOLUME_TYPE::F64;
        _grid = new VolumeGridF64(volumeData->width, volumeData->height, volumeData->depth,
                                  volumeData->data->asDoubles().data());
    } break;

    default:
    {
        LOG_ERROR("Undefined volume format");

    } break;
    }
}
#endif

void Volume::_createGrid(const std::string& hdrFilePath)
{
    // Get the volume meta data from the header file
    const auto volumeData = readHeaderFile(hdrFilePath);

    // Obtain the .img or .bin file path from the prefix and the file type
    std::string parentPath = std::filesystem::path(hdrFilePath).parent_path().string();
    std::string rawFileName = std::filesystem::path(hdrFilePath).stem();
    std::string rawFilePath = parentPath + "/" + rawFileName;

    switch (volumeData->type)
    {
    case VOLUME_TYPE::BIT:
    {
        // Get the raw file path
        rawFilePath += BINARY_EXTENSION;
        // std::string data = readRawFile(rawFilePath);

        std::vector<uint8_t> data = readRawFileToByteVector(rawFilePath);
        std::cout << rawFilePath << std::endl;

        _grid = new BitVolumeGrid(volumeData->width, volumeData->height, volumeData->depth,
                                   new BitArray(data.data(), volumeData->width * volumeData->height * volumeData->depth));

        // Update the grid type
        _gridType = VOLUME_TYPE::BIT;

    } break;

    case VOLUME_TYPE::UI8:
    {
        // Get the raw file path
        rawFilePath += RAW_EXTENSION;
        std::string data = readRawFile(rawFilePath);

        _grid = new VolumeGridU8(volumeData->width, volumeData->height, volumeData->depth,
                                 Array::convertStringTo8UIArray(data));

        // Update the grid type
        _gridType = VOLUME_TYPE::UI8;

    } break;

    case VOLUME_TYPE::UI16:
    {
        // Get the raw file path
        rawFilePath += RAW_EXTENSION;
        std::string data = readRawFile(rawFilePath);

        _grid = new VolumeGridU16(volumeData->width, volumeData->height, volumeData->depth,
                                  Array::convertStringTo16UIArray(data));

        // Update the grid type
        _gridType = VOLUME_TYPE::UI16;

    } break;

    case VOLUME_TYPE::UI32:
    {
        // Get the raw file path
        rawFilePath += RAW_EXTENSION;
        std::string data = readRawFile(rawFilePath);

        _grid = new VolumeGridU32(volumeData->width, volumeData->height, volumeData->depth,
                                  Array::convertStringTo32UIArray(data));

        // Update the grid type
        _gridType = VOLUME_TYPE::UI32;
    } break;

    case VOLUME_TYPE::UI64:
    {
        // Get the raw file path
        rawFilePath += RAW_EXTENSION;
        std::string data = readRawFile(rawFilePath);

        _grid = new VolumeGridU64(volumeData->width, volumeData->height, volumeData->depth,
                                  Array::convertStringTo64UIArray(data));

        // Update the grid type
        _gridType = VOLUME_TYPE::UI64;
    } break;

    case VOLUME_TYPE::F32:
    case VOLUME_TYPE::F64:
    default:
    {
        LOG_ERROR("Unimplemented volume format");
    } break;
    }
}

void Volume::_createGrid(const UltraliserVolumeData* volumeData)
{
    switch (volumeData->type)
    {
    case VOLUME_TYPE::BIT:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::BIT;
        _grid = new BitVolumeGrid(volumeData->width, volumeData->height, volumeData->depth,
                                  Ultraliser::Array::convertStringToBitArray(
                                      volumeData->data,
                                      volumeData->width * volumeData->height * volumeData->depth));
    } break;

    case VOLUME_TYPE::UI8:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI8;
        _grid = new VolumeGridU8(volumeData->width, volumeData->height, volumeData->depth,
                                 Ultraliser::Array::convertStringTo8UIArray(volumeData->data));

    } break;

    case VOLUME_TYPE::UI16:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI16;
        _grid = new VolumeGridU16(volumeData->width, volumeData->height, volumeData->depth,
                                 Ultraliser::Array::convertStringTo16UIArray(volumeData->data));
    } break;

    case VOLUME_TYPE::UI32:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI32;
        _grid = new VolumeGridU8(volumeData->width, volumeData->height, volumeData->depth,
                                 Ultraliser::Array::convertStringTo32UIArray(volumeData->data));

    } break;

    case VOLUME_TYPE::UI64:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::UI64;
        _grid = new VolumeGridU8(volumeData->width, volumeData->height, volumeData->depth,
                                 Ultraliser::Array::convertStringTo64UIArray(volumeData->data));

    } break;

    case VOLUME_TYPE::F32:
    {
        // Create the grid and update the type
        _gridType = VOLUME_TYPE::F32;
    } break;

    case VOLUME_TYPE::F64:
    {

    } break;

    default:
    {
        LOG_ERROR("Undefined volume format");

    } break;
    }
}

void Volume::_allocateGrid(const size_t &width, const size_t &height, const size_t &depth)
{
    // Create the grid
    switch (_gridType)
    {
    case VOLUME_TYPE::BIT:
    {
        _grid = new BitVolumeGrid(width, height, depth);
    } break;

    case VOLUME_TYPE::UI8:
    {
        _grid = new VolumeGridU8(width, height, depth);
    } break;

    case VOLUME_TYPE::UI16:
    {
        _grid = new VolumeGridU16(width, height, depth);
    } break;

    case VOLUME_TYPE::UI32:
    {
        _grid = new VolumeGridU32(width, height, depth);
    } break;

    case VOLUME_TYPE::UI64:
    {
        _grid = new VolumeGridU64(width, height, depth);
    } break;

    case VOLUME_TYPE::F32:
    {
        _grid = new VolumeGridF32(width, height, depth);
    } break;

    case VOLUME_TYPE::F64:
    {
        _grid = new VolumeGridF64(width, height, depth);
    }  break;

    default:
        break;
    }
}

void Volume::_createGrid(void)
{
    LOG_TITLE("Creating Volume Grid");

    // Compute the bounding box size of the given mesh
    Vector3f boundingBoxSize = (_pMax - _pMin);

    // Update the bounding box based on the volume expansion ratio
    if (_expansionRatio > 0.f)
    {
        _pMin -= _expansionRatio * boundingBoxSize;
        _pMax += _expansionRatio * boundingBoxSize;
        boundingBoxSize = (_pMax - _pMin);
    }

    // Find the largest dimension of the mesh model to be able to create a scaled grid.
    _largestDimensionIdx = getLargestDimension(boundingBoxSize);

    // Compute the voxel size
    _voxelSize = boundingBoxSize[_largestDimensionIdx] / (1.f * _baseResolution);

    // Compute the volume dimensions based on the current voxel size
    auto width = F2UI64(std::round(boundingBoxSize[0] / _voxelSize));
    auto height = F2UI64(std::round(boundingBoxSize[1] / _voxelSize));
    auto depth = F2UI64(std::round(boundingBoxSize[2] / _voxelSize));

    LOG_SUCCESS("Volume Dimenions [ %d x %d x %d ] : [ %f x %f x %f ]",
                width, height, depth,
                boundingBoxSize[0], boundingBoxSize[1], boundingBoxSize[2]);

    // Allocating the grid
    _allocateGrid(width, height, depth);
}

void Volume::surfaceVoxelization(Mesh* mesh,
                                 const bool& verbose,
                                 const bool parallel,
                                 const float& sideRatio)
{
    if (verbose) LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    if (verbose)
        LOG_STATUS("Creating Volume Shell [%zu x %zu x %zu]",
                   _grid->getWidth(), _grid->getHeight(), _grid->getDepth());
    if (parallel)
        _rasterizeParallel(mesh, _grid, sideRatio);
    else
        _rasterize(mesh , _grid, sideRatio, verbose);
    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    if (verbose)
    {
        LOG_STATUS_IMPORTANT("Rasterization Stats.");
        LOG_STATS(_surfaceVoxelizationTime);
    }
}

void Volume::surfaceVoxelizationReion(Mesh* mesh,
                                      const Vector3f& pMinRegion,
                                      const Vector3f& pMaxRegion,
                                      const bool& verbose)
{
    if (verbose) LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    if (verbose)
        LOG_STATUS("Creating Volume Shell [%zu x %zu x %zu]",
                   _grid->getWidth(), _grid->getHeight(), _grid->getDepth());

    _rasterizeRegion(mesh , _grid, pMinRegion, pMaxRegion, verbose);
    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    if (verbose)
    {
        LOG_STATUS_IMPORTANT("Rasterization Stats.");
        LOG_STATS(_surfaceVoxelizationTime);
    }
}

void Volume::surfaceVoxelizeNeuronMorphology(
    NeuronMorphology* neuronMorphology, const std::string& packingAlgorithm)
{
    LOG_TITLE("Neuron Surface Voxelization");

    // Start the timer
    TIMER_SET;

    // Get all the sections of the vascular morphology
    Sections sections = neuronMorphology->getSections();

    LOG_STATUS("Creating Volume Shell from Sections");
    LOOP_STARTS("Rasterization");
    PROGRESS_SET;

    // Construct the soma geometry
    auto mesh = new Mesh(neuronMorphology);
    _rasterize(mesh, _grid);

    if (packingAlgorithm == POLYLINE_SPHERE_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); ++i)
        {
            auto section = sections[i];
            auto samples = section->getSamples();

            // If the section contains more than a single sample
            if (samples.size() == 1)
            {
                // Rasterize a single sample
                _rasterize(samples[0], _grid);
            }
            else if (samples.size() > 1)
            {
                // Rasterize a polyline representing the section samples
                auto mesh = new Mesh(samples);
                _rasterize(mesh, _grid);

                // Rasterize the first and last samples as spheres to fill any gaps
                _rasterize(samples.front(), _grid);
                _rasterize(samples.back(), _grid);
            }
            else
            {
                continue;
            }

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else if (packingAlgorithm == POLYLINE_PACKING)
    {
        // Rasterize the first section's beginnigs to fill the gap between the soma and the first
        // sections
        for (auto section: neuronMorphology->getFirstSections())
        {
            auto sample = section->getSamples()[0];
#ifdef ULTRALISER_USE_EIGEN3 
            _rasterize(sample, _grid);
#else
            Samples samples;
            samples.push_back(sample);
            samples.push_back(new Sample(neuronMorphology->getSomaCenter(), sample->getRadius(), 0));
            auto mesh = new Mesh(samples);
            _rasterize(mesh, _grid);
#endif
        }

        Paths paths;
        for (size_t i = 0; i < sections.size(); i++)
        {
            // Construct the paths
            auto section = sections[i];
            Paths sectionPaths = neuronMorphology->getConnectedPathsFromParentsToChildren(section);
            paths.insert(paths.end(), sectionPaths.begin(), sectionPaths.end());
        }

        // Construct the neurites geometry
        OMP_PARALLEL_FOR
                for (size_t i = 0; i < paths.size(); i++)
        {
            auto samples = paths[i];
            auto mesh = new Mesh(samples);
            _rasterize(mesh, _grid);

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, paths.size());
            PROGRESS_UPDATE;
        }
    }
    else if (packingAlgorithm == SDF_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); ++i)
        {
            auto section = sections[i];
            auto samples = section->getSamples();

            // Rasterize each segment of the section samples
            for (size_t i = 0; i < samples.size() - 1; ++i)
                _rasterize(samples[i], samples[i + 1], _grid);

            // Rasterize the last segment of the section
            auto children = section->getChildrenIndices();
            if (children.size() > 0)
                _rasterize(samples.back(), sections[children[0]]->getSamples()[0], _grid);

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else
    {
        LOG_ERROR("[%s] is not a correct packing algorithm.");
    }

    LOOP_DONE;

    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATUS_IMPORTANT("Rasterization Stats.");
    LOG_STATS(_surfaceVoxelizationTime);
}


void Volume::surfaceVoxelizeVasculatureMorphology(
        VasculatureMorphology* vasculatureMorphology,
        const std::string &packingAlgorithm)
{
    LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    // Get all the sections of the vascular morphology
    Sections sections = vasculatureMorphology->getSections();

    LOG_STATUS("Creating Volume Shell from Sections");
    LOOP_STARTS("Rasterization");
    PROGRESS_SET;
    if (packingAlgorithm == POLYLINE_SPHERE_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); ++i)
        {
            auto section = sections[i];
            auto samples = section->getSamples();

           // Rasterize a polyline representing the section samples
            auto mesh = new Mesh(samples);
            _rasterize(mesh, _grid);
            mesh->~Mesh();

            // Rasterize the first and last samples as spheres to fill any gaps
            _rasterize(samples.front(), _grid);
            _rasterize(samples.back(), _grid);

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else if (packingAlgorithm == POLYLINE_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); i++)
        {
            // Construct the paths
            Paths paths = vasculatureMorphology->getConnectedPathsFromParentsToChildren(sections[i]);
            for (size_t j = 0; j < paths.size(); ++j)
            {
                auto mesh = new Mesh(paths[j]);
                _rasterize(mesh , _grid);
                mesh->~Mesh();
            }

            paths.clear();
            paths.shrink_to_fit();

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else if (packingAlgorithm == SDF_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); ++i)
        {
            auto section = sections[i];
            auto samples = section->getSamples();

            // Rasterize each segment of the section samples
            for (uint32_t i = 0; i < samples.size() - 1; ++i)
                _rasterize(samples[i], samples[i + 1], _grid);

            // Rasterize the last segment of the section
            auto children = section->getChildrenIndices();
            if (children.size() > 0)
                _rasterize(samples.back(), sections[children[0]]->getSamples()[0], _grid);
            children.clear();
            children.shrink_to_fit();

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else
    {
        LOG_ERROR("[%s] is not a correct packing algorithm.");
    }
    LOOP_DONE;

    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATUS_IMPORTANT("Rasterization Stats.");
    LOG_STATS(_surfaceVoxelizationTime);
}

void Volume::surfaceVoxelizeAstrocyteMorphology(
        const AstrocyteMorphology* astrocyteMorphology, float threshold, const std::string &packingAlgorithm)
{
    LOG_TITLE("Astrocyte Surface Voxelization");

    // Start the timer
    TIMER_SET;

    // Get all the sections of the vascular morphology
    Sections sections = astrocyteMorphology->getSections();

    LOG_STATUS("Creating Volume Shell from Sections");
    LOOP_STARTS("Rasterization");
    PROGRESS_SET;

    // Use the morphology skeleton to construct a somatic mesh for the astrocyte
    Mesh* somaMesh = new Mesh(astrocyteMorphology);

    // Rasterize the somatic mesh into the volume grid
    _rasterize(somaMesh, _grid);

    if (packingAlgorithm == POLYLINE_SPHERE_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); ++i)
        {
            auto section = sections[i];
            auto samples = section->getSamples();

            // If the section contains more than a single sample
            if (samples.size() == 1)
            {
                // Rasterize a single sample
                _rasterize(samples[0], _grid);
            }
            else if (samples.size() > 1)
            {
                // Rasterize a polyline representing the section samples
                auto mesh = new Mesh(samples);
                _rasterize(mesh, _grid);

                // Rasterize the first and last samples as spheres to fill any gaps
                _rasterize(samples.front(), _grid);
                _rasterize(samples.back(), _grid);
            }
            else
            {
                continue;
            }

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else if (packingAlgorithm == POLYLINE_PACKING)
    {
        // Construct the paths representing the astrocytic processes
        Paths paths;
        for (size_t i = 0; i < sections.size(); i++)
        {
            // For every section, construct a list of paths connecting its parents and children
            Paths sectionPaths =
                    astrocyteMorphology->getConnectedPathsFromParentsToChildren(sections[i]);

            // Add all the paths to the list that will be rasterized
            paths.insert(paths.end(), sectionPaths.begin(), sectionPaths.end());
        }

        // Get reference to the first sections, to handle the connectivity with the soma.
        auto firstSections = astrocyteMorphology->getFirstSections();

        // Construct a straight sections between the somatic center and the neurites
        for (size_t i = 0; i < firstSections.size(); i++)
        {
            // Get the samples of the section
            Samples samples = firstSections[i]->getSamples();

            // Construct a new sample at the center of the soma
            Vector3f newSamplePos = astrocyteMorphology->getSomaCenter();

            // Append the path
            samples.insert(samples.begin(), new Sample(newSamplePos, samples[0]->getRadius(), i));

            // Add the new section to the paths
            paths.push_back(samples);
         }

        // Rasterize all the paths
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < paths.size(); i++)
        {
            // Create a proxy-mesh representing the path
            auto pathProxyMesh = new Mesh(paths[i]);

            // Rasterize the proxy-mesh in the volume grid
            _rasterize(pathProxyMesh, _grid);

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, paths.size());
            PROGRESS_UPDATE;
        }
    }
    else if (packingAlgorithm == SDF_PACKING)
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < sections.size(); ++i)
        {
            auto section = sections[i];
            auto samples = section->getSamples();

            // Rasterize each segment of the section samples
            for (size_t i = 0; i < samples.size() - 1; ++i)
                _rasterize(samples[i], samples[i + 1], _grid);

            // Rasterize the last segment of the section
            auto children = section->getChildrenIndices();
            if (children.size() > 0)
                _rasterize(samples.back(), sections[children[0]]->getSamples()[0], _grid);

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, sections.size());
            PROGRESS_UPDATE;
        }
    }
    else
    {
        LOG_ERROR("[%s] is not a correct packing algorithm.");
    }
    LOOP_DONE;

    // Rasterize the endfeet
    PROGRESS_RESET;
    EndfeetPatches endfeetPatches = astrocyteMorphology->getEndfeetPatches();
    for (size_t j = 0; j < endfeetPatches.size(); ++j)
    {
        EndfootPatches efPatches = endfeetPatches[j];

        // Rasterize sample triangles
        for (size_t i = 0; i < efPatches.size(); ++i)
        {
            // Get the triangle samples
            const Sample* sample0 = efPatches[i]->sample0;
            const Sample* sample1 = efPatches[i]->sample1;
            const Sample* sample2 = efPatches[i]->sample2;

            // Get the vertices positions
            const Vector3f pos0 = sample0->getPosition();
            const Vector3f pos1 = sample1->getPosition();
            const Vector3f pos2 = sample2->getPosition();

            // Get the vertices radii
            const float radius0 = sample0->getRadius();
            const float radius1 = sample1->getRadius();
            const float radius2 = sample2->getRadius();

            // Compute maximum distance based on minimum radius and a threshold
            float minRadius = radius0;
            minRadius = std::min< float >(minRadius, radius1);
            minRadius = std::min< float >(minRadius, radius2);
            float maxSeparation = minRadius * threshold;

            // Compute number of divisions for opposite edge to the sample0
            float edgeDistance = (pos2 - pos1).abs();
            size_t edgeNumDivisions = std::ceil(edgeDistance / maxSeparation);
            float edgeAlphaIncrement = 1.0f / edgeNumDivisions;

            for (size_t j = 0; j < edgeNumDivisions + 1; ++j)
            {
                // Compute opposite edge samples
                float edgeAlpha = j * edgeAlphaIncrement;
                Vector3f edgePos = (1.f - edgeAlpha) * pos1 + edgeAlpha * pos2;
                float edgeRadius = (1.f - edgeAlpha) * radius1 + edgeAlpha * radius2;

                // Compute interpolated samples for the edges betweem the sample0 and the
                // opposite edge samples computed
                float distance = (edgePos - pos0).abs();
                size_t numDivisions = std::ceil(distance / maxSeparation);
                float alphaIncrement = 1.0f / numDivisions;

                // Rasterize sample
                for (size_t k = 0; k < numDivisions + 1; ++k)
                {
                    float alpha = k * alphaIncrement;
                    Vector3f position = (1.0f - alpha) * pos0 + alpha * edgePos;
                    float radius = (1.0f - alpha) * radius0 + alpha * edgeRadius;
                    Sample sample(position, radius, 0);
                    _rasterize(&sample, _grid);
                }
            }
        }
    }

    // Timer
    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATUS_IMPORTANT("Rasterization Stats.");
    LOG_STATS(_surfaceVoxelizationTime);
}

void Volume::surfaceVoxelization(AdvancedMesh *mesh)
{
    LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Creating Volume Shell");
    _rasterize(mesh , _grid);

    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATS(_surfaceVoxelizationTime);
}

void Volume::surfaceVoxelization(const std::string &inputDirectory,
                                 const std::vector< std::string>& meshFiles)
{
    LOG_TITLE("Surface Voxelization");
    TIMER_SET;

    LOG_STATUS("Creating Volume Shell [%d x %d x %d]",
               _grid->getWidth(), _grid->getHeight(), _grid->getDepth());
    size_t processedMeshCount = 0;
    LOOP_STARTS("Rasterization");
    PROGRESS_SET;
    for( size_t iMesh = 0; iMesh < meshFiles.size(); iMesh++ )
    {
        // Create and load the mesh from the file
        std::string meshName = meshFiles[ iMesh ];
        std::string meshFile = meshName;
        if( inputDirectory != EMPTY )
            meshFile = inputDirectory + "/" + meshFile;

        if (File::exists(meshFile))
        {   
            // Input mesh
            auto mesh = new Mesh(meshFile, false);

            // Surface voxelization
            surfaceVoxelization(mesh, false, true);
            processedMeshCount++;

            // Free the mesh
            delete mesh;
        }

        // Update the progress bar
        // LOOP_PROGRESS(PROGRESS, meshFiles.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;

    // Statistics
    _surfaceVoxelizationTime = GET_TIME_SECONDS;
    LOG_STATS(_surfaceVoxelizationTime);
    LOG_DETAIL("[%zu/%zu] Meshes were Voxelized with Surface Voxelization",
               processedMeshCount, meshFiles.size());
}



std::vector< Vec3ui_64 > Volume::verifyBorderVoxels(Mesh* mesh,
                                const float& sideRatio, const bool& verbose)
{

    std::vector< Vec3ui_64 > candidates;

    if (verbose) LOOP_STARTS("Rasterization");
    size_t progress = 0;

    size_t allVoxels = 0, toBeRemoved = 0;

    for (size_t triangleIdx = 0; triangleIdx < mesh->getNumberTriangles(); ++triangleIdx)
    {
        ++progress;
        if (verbose) LOOP_PROGRESS(progress, mesh->getNumberTriangles());

        // Get the pMin and pMax of the triangle within the grid
        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getBoundingBox(mesh, triangleIdx, pMinTriangle, pMaxTriangle);


        for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ++ix)
        {
            for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; ++iy)
            {
                for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; ++iz)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));

                    auto fullTest = _testTriangleCubeIntersection(mesh, triangleIdx, gi);

                    if (fullTest)
                    {
                        auto sideTest = _testTriangleCubeIntersection(mesh, triangleIdx, gi, sideRatio);

                        if (!sideTest)
                        {
                            _grid->clearVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                        }
                        else
                            _grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));

                    }

                }
            }
        }
    }
    if (verbose) LOOP_DONE;

    std::cout << computeNumberNonZeroVoxels() << "\n";

    return candidates;
}

void Volume::_rasterize(Mesh* mesh, VolumeGrid* grid, const float& sideRatio, const bool& verbose)
{
    if (verbose) LOOP_STARTS("Rasterization");
    size_t progress = 0;
    
    for (size_t triangleIdx = 0; triangleIdx < mesh->getNumberTriangles(); ++triangleIdx)
    {
        ++progress;
        if (verbose) LOOP_PROGRESS(progress, mesh->getNumberTriangles());

        // Get the pMin and pMax of the triangle within the grid
        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getBoundingBox(mesh, triangleIdx, pMinTriangle, pMaxTriangle);

        for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ++ix)
        {
            for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; ++iy)
            {
                for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; ++iz)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                    if (_testTriangleCubeIntersection(mesh, triangleIdx, gi, sideRatio))
                        grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                }
            }
        }
    }
    if (verbose) LOOP_DONE;
}


void Volume::_rasterizeRegion(Mesh* mesh, VolumeGrid* grid,
                              const Vector3f& pMinRegion,
                              const Vector3f& pMaxRegion,
                              const bool& verbose)
{
    if (verbose) LOOP_STARTS("Rasterization");
    size_t progress = 0;

    for (size_t triangleIdx = 0; triangleIdx < mesh->getNumberTriangles(); ++triangleIdx)
    {
        ++progress;
        if (verbose) LOOP_PROGRESS(progress, mesh->getNumberTriangles());

        // Get the pMin and pMax of the triangle within the grid
        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getBoundingBox(mesh, triangleIdx, pMinTriangle, pMaxTriangle);


        // Get the indices of the grid that correspond to the bounding box
        GridIndex vMin, vMax;
        _vec2grid(pMinRegion, vMin);
        _vec2grid(pMaxRegion, vMax);


        int64_t tMin[3];
        int64_t tMax[3];

        tMin[0] = vMin[0]; tMin[1] = vMin[1]; tMin[2] = vMin[2];
        tMax[0] = vMax[0]; tMax[1] = vMax[1]; tMax[2] = vMax[2];

        tMin[0] = std::max(int64_t(0), tMin[0]);
        tMax[0] = std::min(int64_t(_grid->getWidth() - 1) , tMax[0]);

        tMin[1] = std::max(int64_t(0), tMin[1]);
        tMax[1] = std::min(int64_t(_grid->getHeight() - 1) , tMax[1]);

        tMin[2] = std::max(int64_t(0), tMin[2]);
        tMax[2] = std::min(int64_t(_grid->getDepth() - 1) , tMax[2]);


        // std::cout << pMinTriangle[0] << " " << pMinTriangle[1] << " " << pMinTriangle[2] << "\n";
        if (pMinTriangle[0] >= tMin[0] &&
            pMinTriangle[1] >= tMin[1] &&
            pMinTriangle[2] >= tMin[2] &&
            pMaxTriangle[0] <= tMax[0] &&
            pMaxTriangle[1] <= tMax[1] &&
            pMaxTriangle[2] <= tMax[2])
        {
            for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ++ix)
            {
                for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; ++iy)
                {
                    for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; ++iz)
                    {
                        GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                        if (_testTriangleCubeIntersection(mesh, triangleIdx, gi))
                            grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                    }
                }
            }
        }
    }
    if (verbose) LOOP_DONE;
}

void Volume::_rasterize(Sample* sample, VolumeGrid* grid)
{
    // Get the pMin and pMax of the sphere within the grid
    int64_t pMinTriangle[3], pMaxTriangle[3];
    _getBoundingBox(sample, pMinTriangle, pMaxTriangle);

    for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ++ix)
    {
        for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; ++iy)
        {
            for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; ++iz)
            {
                GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                if (_testSampleCubeIntersection(sample, gi))
                    grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
            }
        }
    }
}

void Volume::_rasterize(Sample* sample0, Sample* sample1, VolumeGrid* grid, float stepAlpha)
{
    Vector3f position0 = sample0->getPosition();
    Vector3f position1 = sample1->getPosition();
    float radius0 = sample0->getRadius();
    float radius1 = sample1->getRadius();

    // Compute the number of interpolation steps, and the position and radius increment    
    float minRadius = std::min(radius0, radius1) * stepAlpha;
    float distance = (position1 - position0).abs();
    uint32_t numSteps = std::ceil( distance/ minRadius);
    Vector3f positionIncrement = (position1 - position0) / numSteps;
    float radiusIncrement =  (radius1 - radius0) / numSteps;

    for (uint32_t i = 0; i <= numSteps; ++i)
    {
        // Compute the interpolated sample
        Sample sample(position0 + positionIncrement * i,
                      radius0 + radiusIncrement * i, 0);
        _rasterize(&sample, grid);
    }
}

void Volume::_rasterizeParallel(Mesh* mesh, VolumeGrid* grid, const float& sideRatio)
{
    // Start the timer
    TIMER_SET;

    LOOP_STARTS("Rasterization");
    PROGRESS_SET;
    #pragma omp parallel for schedule(dynamic)
    for (size_t tIdx = 0; tIdx < mesh->getNumberTriangles(); tIdx++)
    {
        // Get the pMin and pMax of the triangle within the grid
        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getBoundingBox(mesh, tIdx, pMinTriangle, pMaxTriangle);

        for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ix++)
        {
            for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; iy++)
            {
                for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; iz++)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                    if (_testTriangleCubeIntersection(mesh, tIdx, gi, sideRatio))
                        grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                }
            }
        }

        // Update the progress bar
         LOOP_PROGRESS(PROGRESS, mesh->getNumberTriangles());
         PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::_rasterize(AdvancedMesh* mesh, VolumeGrid* grid)
{
    // Get a an array of triangles
    AdvancedTriangle** triangles = (AdvancedTriangle **) mesh->_triangles.toArray();
    int triangleCount = mesh->_triangles.numberElements();

    LOOP_STARTS("Rasterization");
    size_t progress = 0;
    for (int i = 0; i <triangleCount; ++i)
    {
        AdvancedTriangle triangle = *triangles[i];

        ++progress;
        LOOP_PROGRESS(progress, triangleCount);

        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getTriangleBoundingBox(triangle, pMinTriangle, pMaxTriangle);

        for (int64_t ix = pMinTriangle[0];
             ix <= (pMaxTriangle[0]); ix++)
        {
            for (int64_t iy = pMinTriangle[1];
                 iy <= (pMaxTriangle[1]); iy++)
            {
                for (int64_t iz = pMinTriangle[2];
                     iz <= (pMaxTriangle[2]); iz++)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                    if (_testTriangleGridIntersection(triangle, gi))
                        grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                }
            }
        }
    }
    LOOP_DONE;
}

void Volume::solidVoxelization(const SOLID_VOXELIZATION_AXIS& axis)
{
    LOG_TITLE("Solid Voxelization");

    // The 2D flood filling is only supported for the solid voxelization
    LOG_STATUS("Filling Volume");
    _floodFill2D(axis);

    LOG_STATUS_IMPORTANT("Solid Voxelization Stats.");
    LOG_STATS(_solidVoxelizationTime);
}

void Volume::_floodFill2D(const SOLID_VOXELIZATION_AXIS &axis)
{
    // Start the timer
    TIMER_SET;

    switch (axis)
    {
    case X: _floodFillAlongAxis(_grid, SOLID_VOXELIZATION_AXIS::X);
        break;

    case Y: _floodFillAlongAxis(_grid, SOLID_VOXELIZATION_AXIS::Y);
        break;

    case Z: _floodFillAlongAxis(_grid, SOLID_VOXELIZATION_AXIS::Z);
        break;

    case XYZ:_floodFillAlongXYZ(_grid);
        break;
    }

    // Save the solid voxelization time
    _solidVoxelizationTime = GET_TIME_SECONDS;
}

void Volume::_floodFillAlongAxis(VolumeGrid* grid, const SOLID_VOXELIZATION_AXIS &axis)
{
    /// Disable buffering
    setbuf(stdout, nullptr);

    // Start the timer
    TIMER_SET;

    // The dimension with which the flood filling will happen
    int64_t dimension;

    // Flood-filling string
    std::string floodFillingString;

    // Flood-filling axis
    AXIS floodFillingAxis;

    switch (axis)
    {
    case SOLID_VOXELIZATION_AXIS::X:
    {
        dimension = getWidth();
        floodFillingAxis = AXIS::X;
        floodFillingString = "2D Slice Flood-filling (X-axis)";
    } break;

    case SOLID_VOXELIZATION_AXIS::Y:
    {
        dimension = getHeight();
        floodFillingAxis = AXIS::Y;
        floodFillingString = "2D Slice Flood-filling (Y-axis)";
    } break;

    case SOLID_VOXELIZATION_AXIS::Z:
    {
        dimension = getDepth();
        floodFillingAxis = AXIS::Z;
        floodFillingString = "2D Slice Flood-filling (Z-axis)";
    } break;

    // XYZ voxelization will be handled
    case SOLID_VOXELIZATION_AXIS::XYZ:
        break;
    }

    LOOP_STARTS(floodFillingString.c_str());
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = 0 ; i < dimension; ++i)
    {
        grid->floodFillSliceAlongAxis(i, floodFillingAxis);

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, dimension);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::_floodFillAlongXYZ(VolumeGrid *grid)
{
    // Volume grids per axis
    VolumeGrid *xGrid, *yGrid, *zGrid;

    // Create the grid
    switch (_gridType)
    {
    case VOLUME_TYPE::BIT:
    {
        xGrid = new BitVolumeGrid(static_cast< BitVolumeGrid* >(grid));
        yGrid = new BitVolumeGrid(static_cast< BitVolumeGrid* >(grid));
        zGrid = new BitVolumeGrid(static_cast< BitVolumeGrid* >(grid));

    } break;

    case VOLUME_TYPE::UI8:
    {
        xGrid = new UnsignedVolumeGrid<uint8_t>(static_cast< UnsignedVolumeGrid<uint8_t>* >(grid));
        yGrid = new UnsignedVolumeGrid<uint8_t>(static_cast< UnsignedVolumeGrid<uint8_t>* >(grid));
        zGrid = new UnsignedVolumeGrid<uint8_t>(static_cast< UnsignedVolumeGrid<uint8_t>* >(grid));
    } break;

    case VOLUME_TYPE::UI16:
    {
        xGrid = new UnsignedVolumeGrid<uint16_t>(static_cast< UnsignedVolumeGrid<uint16_t>* >(grid));
        yGrid = new UnsignedVolumeGrid<uint16_t>(static_cast< UnsignedVolumeGrid<uint16_t>* >(grid));
        zGrid = new UnsignedVolumeGrid<uint16_t>(static_cast< UnsignedVolumeGrid<uint16_t>* >(grid));
    } break;

    case VOLUME_TYPE::UI32:
    {
        xGrid = new UnsignedVolumeGrid<uint32_t>(static_cast< UnsignedVolumeGrid<uint32_t>* >(grid));
        yGrid = new UnsignedVolumeGrid<uint32_t>(static_cast< UnsignedVolumeGrid<uint32_t>* >(grid));
        zGrid = new UnsignedVolumeGrid<uint32_t>(static_cast< UnsignedVolumeGrid<uint32_t>* >(grid));
    } break;

    case VOLUME_TYPE::UI64:
    {
        xGrid = new UnsignedVolumeGrid<uint64_t>(static_cast< UnsignedVolumeGrid<uint64_t>* >(grid));
        yGrid = new UnsignedVolumeGrid<uint64_t>(static_cast< UnsignedVolumeGrid<uint64_t>* >(grid));
        zGrid = new UnsignedVolumeGrid<uint64_t>(static_cast< UnsignedVolumeGrid<uint64_t>* >(grid));
    } break;

    case VOLUME_TYPE::F32:
    case VOLUME_TYPE::F64:
    {
        LOG_ERROR("_floodFillAlongXYZ CANNOT be applied to Float Grids!");
    } break;
    }

    // Flood fill along the three axes
    _floodFillAlongAxis(xGrid, SOLID_VOXELIZATION_AXIS::X);
    _floodFillAlongAxis(yGrid, SOLID_VOXELIZATION_AXIS::Y);
    _floodFillAlongAxis(zGrid, SOLID_VOXELIZATION_AXIS::Z);

    // Blend the three grids using AND operation and store the final result in the xGrid
    xGrid->andWithAnotherGrid(yGrid);
    xGrid->andWithAnotherGrid(zGrid);

    // Blend the xGrid with the default _grid
    _grid->orWithAnotherGrid(xGrid);

    // Release the auxiliary grids
    delete xGrid;
    delete yGrid;
    delete zGrid;
}

void Volume::getVoxelBoundingBox(const int64_t& x, const int64_t& y, const int64_t& z,
                                 Vector3f& pMin, Vector3f& pMax) const
{
    // pMin
    pMin.x() = _pMin[0] + (x * _voxelSize);
    pMin.y() = _pMin[1] + (y * _voxelSize);
    pMin.z() = _pMin[2] + (z * _voxelSize);

    // Just add the voxel size along the three-dimensions
    pMax = pMin + Vector3f(_voxelSize);
}

void Volume::getVolumeBoundingBox(Vector3f& pMin, Vector3f& pMax) const
{
    pMin.x() = _pMin.x();
    pMin.y() = _pMin.y();
    pMin.z() = _pMin.z();

    pMax.x() = _pMax.x();
    pMax.y() = _pMax.y();
    pMax.z() = _pMax.z();
}

int Volume::_triangleCubeSign(Mesh *mesh,
                              int tIdx, const GridIndex & gi)
{
    Vector3f boxcenter((0.5f + gi[0]) * _voxelSize + _pMin[0],
                       (0.5f + gi[1]) * _voxelSize + _pMin[1],
                       (0.5f + gi[2]) * _voxelSize + _pMin[2]);

    Vector3f tv[3];
    for (size_t ii = 0; ii < 3; ++ii)
        tv[ii] = mesh->getVertices()[mesh->getTriangles()[tIdx][ii]];

    Vector3f e1 = tv[1] - tv[0];
    Vector3f e2 = tv[2] - tv[0];
    Vector3f n = Vector3f::cross(e1, e2).normalized();
    Vector3f d = boxcenter - tv[0];
    float dotp = Vector3f::dot(n, d);
    if (dotp > 0)
        return 2;

    // Too far away
    if (dotp < -0.9f * _voxelSize)
        return 3;

    n.normalize();
    d = boxcenter - (Vector3f::dot(n, d)) * n;
    Vector3f n0 = Vector3f::cross(tv[1] - d, tv[2] - d);

    // n0.normalize();
    float thresh = 1.0;

    if (Vector3f::dot(n0, n) < -thresh)
        return 1;

    Vector3f n1 = Vector3f::cross(tv[2] - d, tv[0] - d);
    // n1.normalize();
    if (Vector3f::dot(n1, n) < -thresh)
        return 1;

    Vector3f n2 = Vector3f::cross(tv[0] - d, tv[1] - d);
    // n2.normalize();
    if (Vector3f::dot(n2, n)< -thresh)
        return 1;

    return -1;
}

bool Volume::_testTriangleCubeIntersection(Mesh* mesh, size_t triangleIdx,
                                           const GridIndex& voxel, const float &sideRatio)
{
    // Get the origin of the voxel
    double  voxelOrigin[3];
    voxelOrigin[0] = _pMin[0] + (voxel[0] * _voxelSize);
    voxelOrigin[1] = _pMin[1] + (voxel[1] * _voxelSize);
    voxelOrigin[2] = _pMin[2] + (voxel[2] * _voxelSize);

    // Voxel half size
    double voxelHalfSize[3];
    voxelHalfSize[0] = _voxelSize * 0.5;
    voxelHalfSize[1] = _voxelSize * 0.5;
    voxelHalfSize[2] = _voxelSize * 0.5;

    // Voxel center
    double voxelCenter[3];
    voxelCenter[0] = voxelOrigin[0] + voxelHalfSize[0];
    voxelCenter[1] = voxelOrigin[1] + voxelHalfSize[1];
    voxelCenter[2] = voxelOrigin[2] + voxelHalfSize[2];

//    auto v1 = mesh->getVertices()[mesh->getTriangles()[triangleIdx][0]];
//    auto v2 = mesh->getVertices()[mesh->getTriangles()[triangleIdx][1]];
//    auto v3 = mesh->getVertices()[mesh->getTriangles()[triangleIdx][2]];

    // Triangle vertices
    double triangle[3][3];

//    triangle[0][0] = v1.x();
//    triangle[0][1] = v1.y();
//    triangle[0][2] = v1.z();

//    triangle[1][0] = v2.x();
//    triangle[1][1] = v2.y();
//    triangle[1][2] = v2.z();

//    triangle[2][0] = v3.x();
//    triangle[2][1] = v3.y();
//    triangle[2][2] = v3.z();

    // For each vertex in the triangle
    for (size_t i = 0; i < 3; ++i)
    {
        // For each coordinate of the vertex
        for (size_t j = 0; j < 3; ++j)
        {
            // Load all the verticies of the selected triangle in _triangle_
            triangle[i][j] = mesh->getVertices()[mesh->getTriangles()[triangleIdx][i]][j];
        }
    }

    voxelHalfSize[0] *= sideRatio;
    voxelHalfSize[1] *= sideRatio;
    voxelHalfSize[2] *= sideRatio;

    return checkTriangleBoxIntersection(voxelCenter, voxelHalfSize, triangle);
}

bool Volume::_testSampleCubeIntersection(Sample* sample, const GridIndex& voxel)
{
    // Voxel half size
    float voxelHalfSize = _voxelSize * 0.5f;

    // Get the center of the voxel
    Vector3f voxelCenter;
    voxelCenter[0] = _pMin[0] + (voxel[0] * _voxelSize) + voxelHalfSize;
    voxelCenter[1] = _pMin[1] + (voxel[1] * _voxelSize) + voxelHalfSize;
    voxelCenter[2] = _pMin[2] + (voxel[2] * _voxelSize) + voxelHalfSize;
 
    Vector3f position = sample->getPosition() - voxelCenter;
    position[0] = abs(position[0]);
    position[1] = abs(position[1]);
    position[2] = abs(position[2]);
    float r2 = pow(sample->getRadius(), 2);
    float minDist = 0.0f;
    
    for (uint8_t i = 0; i < 3; ++i)
    {
        if (position[i] > voxelHalfSize) minDist += pow(position[i] - voxelHalfSize, 2);
    }
    return minDist <= r2;
}

bool Volume::_testTriangleGridIntersection(AdvancedTriangle triangle,
                                           const GridIndex& voxel)
{
    // Get the origin of the voxel
    double  voxelOrigin[3];
    voxelOrigin[0] = _pMin[0] + (voxel[0] * _voxelSize);
    voxelOrigin[1] = _pMin[1] + (voxel[1] * _voxelSize);
    voxelOrigin[2] = _pMin[2] + (voxel[2] * _voxelSize);

    // Voxel half size
    double voxelHalfSize[3];
    voxelHalfSize[0] = _voxelSize * 0.5;
    voxelHalfSize[1] = _voxelSize * 0.5;
    voxelHalfSize[2] = _voxelSize * 0.5;

    // Voxel center
    double voxelCenter[3];
    voxelCenter[0] = voxelOrigin[0] + voxelHalfSize[0];
    voxelCenter[1] = voxelOrigin[1] + voxelHalfSize[1];
    voxelCenter[2] = voxelOrigin[2] + voxelHalfSize[2];

    // Construct an array to represent the data
    double triangleArray[3][3];
    triangleArray[0][0] = triangle.v1()->x;
    triangleArray[0][1] = triangle.v1()->y;
    triangleArray[0][2] = triangle.v1()->z;

    triangleArray[1][0] = triangle.v2()->x;
    triangleArray[1][1] = triangle.v2()->y;
    triangleArray[1][2] = triangle.v2()->z;

    triangleArray[2][0] = triangle.v3()->x;
    triangleArray[2][1] = triangle.v3()->y;
    triangleArray[2][2] = triangle.v3()->z;

    // Test if the triangle and the voxel are intersecting or not
    return checkTriangleBoxIntersection(voxelCenter, voxelHalfSize, triangleArray);
}

size_t Volume::_clampIndex(size_t idx, size_t dimension)
{
    idx = std::max(static_cast< size_t >(0) , idx);
    idx = std::min(idx, static_cast< size_t >(_grid->getDimension(dimension) - 1));
    return idx;
}

void Volume::_vec2grid(const Vector3f& point, GridIndex& gridIndex)
{
    gridIndex[0] = F2I64((point[0] - _pMin[0]) / _voxelSize);
    gridIndex[1] = F2I64((point[1] - _pMin[1]) / _voxelSize);
    gridIndex[2] = F2I64((point[2] - _pMin[2]) / _voxelSize);
}

void Volume::_getTriangleBoundingBox(AdvancedTriangle triangle, int64_t *tMin, int64_t *tMax)
{
    // Find the index of the voxel that intersects the triangle
    GridIndex vIdx;
    Vertex vertex(triangle.v1()->x, triangle.v1()->y, triangle.v1()->z);

    _vec2grid(vertex, vIdx);

    for (int64_t j = 0; j < DIMENSIONS; ++j)
    {
        tMin[j] = (vIdx[I2UI64(j)] - 1);
        tMax[j] = (vIdx[I2UI64(j)]);
    }

    for (size_t j = 1; j < DIMENSIONS; ++j)
    {
        if (j == 1)
        {
            vertex.x() = triangle.v2()->x;
            vertex.y() = triangle.v2()->y;
            vertex.z() = triangle.v2()->z;
        }
        else
        {
            vertex.x() = triangle.v3()->x;
            vertex.y() = triangle.v3()->y;
            vertex.z() = triangle.v3()->z;
        }

        _vec2grid(vertex, vIdx);

        for (size_t k = 0; k < DIMENSIONS; ++k)
        {
            if (vIdx[k] - 1 < (tMin[k]))
                tMin[k] = (vIdx[k] - 1);

            if (vIdx[k] > (tMax[k]))
                tMax[k] = (vIdx[k]);
        }
    }

    for (int32_t ii = 0; ii < DIMENSIONS ; ++ii)
    {
        tMin[ii] = std::max(int64_t(0), tMin[ii] - 1);
        tMax[ii] = std::min(int64_t(_grid->getDimension(ii) - 2), tMax[ii] + 2);
    }
}

void Volume::_getBoundingBox(Mesh* mesh, size_t i, int64_t *tMin, int64_t *tMax)
{
    // The bounding box of the triangle
    Vector3f pMin, pMax;
    mesh->getTriangleBoundingBox(i, pMin, pMax);

    // Get the indices of the grid that correspond to the bounding box
    GridIndex vMin, vMax;
    _vec2grid(pMin, vMin);
    _vec2grid(pMax, vMax);

    tMin[0] = vMin[0]; tMin[1] = vMin[1]; tMin[2] = vMin[2];
    tMax[0] = vMax[0]; tMax[1] = vMax[1]; tMax[2] = vMax[2];

    tMin[0] = std::max(int64_t(0), tMin[0]);
    tMax[0] = std::min(int64_t(_grid->getWidth() - 1) , tMax[0]);

    tMin[1] = std::max(int64_t(0), tMin[1]);
    tMax[1] = std::min(int64_t(_grid->getHeight() - 1) , tMax[1]);

    tMin[2] = std::max(int64_t(0), tMin[2]);
    tMax[2] = std::min(int64_t(_grid->getDepth() - 1) , tMax[2]);
}

void Volume::_getBoundingBox(Sample* sample, int64_t *tMin, int64_t *tMax)
{
    // The bounding box of the sphere
    Vector3f pMin, pMax;
    Vector3f radius(sample->getRadius());
    pMin = sample->getPosition() - radius;
    pMax = sample->getPosition() + radius;

    // Get the indices of the grid that correspond to the bounding box
    GridIndex vMin, vMax;
    _vec2grid(pMin, vMin);
    _vec2grid(pMax, vMax);

    tMin[0] = vMin[0]; tMin[1] = vMin[1]; tMin[2] = vMin[2];
    tMax[0] = vMax[0]; tMax[1] = vMax[1]; tMax[2] = vMax[2];

    tMin[0] = std::max(int64_t(0), tMin[0]);
    tMax[0] = std::min(int64_t(_grid->getWidth() - 1) , tMax[0]);

    tMin[1] = std::max(int64_t(0), tMin[1]);
    tMax[1] = std::min(int64_t(_grid->getHeight() - 1) , tMax[1]);

    tMin[2] = std::max(int64_t(0), tMin[2]);
    tMax[2] = std::min(int64_t(_grid->getDepth() - 1) , tMax[2]);
}

void Volume::_getBoundingBox(Sample* sample0, Sample* sample1, int64_t *tMin, int64_t *tMax)
{
    // The bounding box of the sphere
    Vector3f pMin, pMax;
    Vector3f radius0(sample0->getRadius());
    Vector3f minimum0 = sample0->getPosition() - radius0;
    Vector3f maximum0 = sample0->getPosition() + radius0;
    Vector3f radius1(sample1->getRadius());
    Vector3f minimum1 = sample1->getPosition() - radius1;
    Vector3f maximum1 = sample1->getPosition() + radius1;
    
    pMin[0] = std::min(minimum0[0], minimum1[0]);
    pMin[1] = std::min(minimum0[1], minimum1[1]);
    pMin[2] = std::min(minimum0[2], minimum1[2]);
    pMax[0] = std::max(maximum0[0], maximum1[0]);
    pMax[1] = std::max(maximum0[1], maximum1[1]);
    pMax[2] = std::max(maximum0[2], maximum1[2]);

    // Get the indices of the grid that correspond to the bounding box
    GridIndex vMin, vMax;
    _vec2grid(pMin, vMin);
    _vec2grid(pMax, vMax);

    tMin[0] = vMin[0]; tMin[1] = vMin[1]; tMin[2] = vMin[2];
    tMax[0] = vMax[0]; tMax[1] = vMax[1]; tMax[2] = vMax[2];

    tMin[0] = std::max(int64_t(0), tMin[0]);
    tMax[0] = std::min(int64_t(_grid->getWidth() - 1) , tMax[0]);

    tMin[1] = std::max(int64_t(0), tMin[1]);
    tMax[1] = std::min(int64_t(_grid->getHeight() - 1) , tMax[1]);

    tMin[2] = std::max(int64_t(0), tMin[2]);
    tMax[2] = std::min(int64_t(_grid->getDepth() - 1) , tMax[2]);
}

int64_t Volume::getWidth(void) const
{
    return _grid->getWidth();
}

int64_t Volume::getHeight(void) const
{
    return _grid->getHeight();
}

int64_t Volume::getDepth(void) const
{
    return _grid->getDepth();
}

size_t Volume::getNumberVoxels(void) const
{
    return static_cast< size_t >(_grid->getWidth() * _grid->getHeight() * _grid->getDepth());
}

double Volume::getSurfaceVoxelizationTime(void) const
{
    return _surfaceVoxelizationTime;
}

double Volume::getSolidVoxelizationTime(void) const
{
    return _solidVoxelizationTime;
}

double Volume::getAppendingVolumePassTime(void) const
{
    return _addingVolumePassTime;
}

bool Volume::isFilled(const u_int64_t& index) const
{
    return _grid->isFilled(index);
}

bool Volume::isFilled(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapToIndex(x, y, z, outlier);
    if (outlier)
        return false;
    else
        return isFilled(index);
}

bool Volume::isFilledWithoutBoundCheck(const int64_t &x, const int64_t &y, const int64_t &z) const
{
   return isFilled(mapTo1DIndexWithoutBoundCheck(x, y, z));
}

size_t Volume::mapTo1DIndexWithoutBoundCheck(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    return I2UI64(x + (getWidth() * y) + (getWidth() * getHeight() * z));
}

size_t Volume::mapToIndex(const int64_t &x, const int64_t &y, const int64_t &z, bool& outlier) const
{
    if(x >= getWidth()  || x < 0 || y >= getHeight() || y < 0 || z >= getDepth()  || z < 0)
    {
        outlier = true;
        return 0;
    }
    else
    {
        outlier = false;
        return I2UI64(x + (getWidth() * y) + (getWidth() * getHeight() * z));
    }
}

void Volume::mapToXYZ(const size_t index, size_t& x, size_t& y, size_t& z) const
{
    size_t idx = index;
    z = idx / (getWidth() * getHeight());
    idx -= z * getWidth() * getHeight();

    y = idx / getWidth();
    idx -= y * getWidth();
          x = idx/1;
          return;

    x = index / (getHeight() * getDepth());
    y = (index / getDepth()) % getHeight();
    z = index % getDepth();
}

void Volume::project(const std::string prefix,
                     const bool xy, const bool xz, const bool zy,
                     const bool &projectColorCoded) const
{
    _grid->projectVolume(prefix, xy, xz, zy, projectColorCoded);
}

void Volume::writeVolumes(const std::string &prefix,
                          const bool& bitFormat, const bool &unsignedFormat, const bool &floatFormat,
                          const bool& nrrdFormat, const bool& rawFormat) const
{
    if (bitFormat || unsignedFormat || floatFormat || nrrdFormat || rawFormat)
    {
        // Start the timer
        TIMER_SET;

        LOG_TITLE("Writing Volumes");

        if (bitFormat)
        {
            LOG_SUCCESS("Bit Volume (1-bit per voxel)");
            _grid->writeBitVolume(prefix);
        }

        if (unsignedFormat)
        {
            LOG_SUCCESS("Unsigned Volume");
            _grid->writeUnsignedVolume(prefix);
        }

        if (floatFormat)
        {
            LOG_SUCCESS("Float Volume");
            _grid->writeFloatVolume(prefix);
        }

        if (nrrdFormat)
        {
            LOG_SUCCESS("NRRD Raw Volume in .NRRD file");
            _grid->writeNRRDVolume(prefix);
        }

        if (rawFormat)
        {
            LOG_SUCCESS("Raw Volume (1, 2, 3 or 4-bytes per voxel in .HDR/.IMG files)");
            _grid->writeRAWVolume(prefix);
        }

        // Statictics
        LOG_STATUS_IMPORTANT("Writing Volumes Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void Volume::writeStackXY(const std::string &outputDirectory, const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "xy-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: XY");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t z = 0; z < getDepth(); ++z)
    {
        // Create a slice
        auto slice = std::make_unique<Image>(getWidth(), getHeight());

        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                bool outlier;
                size_t index = mapToIndex(I2I64(i), I2I64(j), I2I64(getDepth() - 1 - z), outlier);

                if (_grid->isFilled(index) && !outlier)
                    slice->setPixelColor(i , j, WHITE);
                else
                    slice->setPixelColor(i , j, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << z;
        slice->writePPM(stream.str());

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, getDepth());
        PROGRESS_UPDATE;
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStackXZ(const std::string &outputDirectory,
                            const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "xz-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: XZ");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t y = 0; y < getHeight(); ++y)
    {
        // Create a slice
        auto slice = std::make_unique<Image>(getWidth(), getHeight());

        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t k = 0; k < getDepth(); k++)
            {
                bool outlier;
                size_t index = mapToIndex(I2I64(i), I2I64(k), I2I64(getHeight() - 1 - y), outlier);

                if (_grid->isFilled(index) && !outlier)
                    slice->setPixelColor(i , k, WHITE);
                else
                    slice->setPixelColor(i , k, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << y;
        slice->writePPM(stream.str());

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, getHeight());
        PROGRESS_UPDATE;
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStackZY(const std::string &outputDirectory,
                            const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "zy-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: ZY");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = 0; i < getWidth(); ++i)
    {
        auto slice = std::make_unique<Image>(getDepth(), getHeight());

        for (int64_t z = 0; z < getDepth(); z++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                bool outlier;
                size_t index = mapToIndex(I2I64(getWidth() - 1 - i), I2I64(j), I2I64(z), outlier);

                if (_grid->isFilled(index) && !outlier)
                    slice->setPixelColor(z , getHeight() - j - 1, WHITE);
                else
                    slice->setPixelColor(z , getHeight() - j - 1, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << i;
        slice->writePPM(stream.str());

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, getWidth());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStacks(const std::string &outputDirectory,
                         const std::string &prefix,
                         const bool& xy,
                         const bool &xz,
                         const bool& zy) const
{
    // Start timer
    TIMER_SET;

    if (xy || xz || zy)
    {
        LOG_TITLE("Writing Stacks");

        if (xy)
            writeStackXY(outputDirectory, prefix);

        if (xz)
            writeStackXZ(outputDirectory, prefix);

        if (zy)
            writeStackZY(outputDirectory, prefix);

        // Statictics
        LOG_STATUS_IMPORTANT("Writing Stacks Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void Volume::exportToMesh(const std::string &prefix,
                          const bool &formatOBJ, const bool &formatPLY,
                          const bool &formatOFF, const bool &formatSTL) const
{
    if (!(formatOBJ || formatPLY || formatOFF || formatSTL))
    {
        LOG_WARNING("Exporto mesh option must be enabled to export this mesh. "
                    "User one of the following: "
                    "[--export-obj, --export-ply, --export-off, --export-stl]");
        return;
    }

    TIMER_SET;
    LOG_TITLE("Constructing Volume Mesh");

    // The generated mesh from the volume
    std::unique_ptr< VolumeMesh > volumeMesh = std::make_unique< VolumeMesh >();

    // Delta value
    const Vector3f delta(1, 1, 1);

    LOOP_STARTS("Searching Filled Voxelss")
    for (int64_t i = 0; i < _grid->getWidth(); ++i)
    {
        LOOP_PROGRESS(i, _grid->getWidth());

        for (int64_t j = 0; j < _grid->getHeight(); ++j)
        {
            for (int64_t k = 0; k < _grid->getDepth(); ++k)
            {
                // Skip empty voxels
                if (_grid->isEmpty(i, j, k))
                    continue;

                Vector3f coordinate(i, j, k);
                Vector3f pMin = _pMin + (_voxelSize * coordinate);
                Vector3f pMax = pMin + Vector3f(_voxelSize);

                // A mesh representing the bounding box of the cube
                VolumeMesh* voxelCube = VolumeMesh::constructVoxelCube(pMin, pMax);

                // Append it to the volume mesh
                volumeMesh->append(voxelCube);

                // Free the voxel cube
                voxelCube->~VolumeMesh();
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOG_STATUS_IMPORTANT("Volume Mesh Construction Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    LOG_TITLE("Exporting Volume Mesh");
    TIMER_RESET;
    const std::string outputPrefix = prefix + VOLUME_MESH_SUFFIX;
    if (formatOBJ)
    {
        exportOBJ(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatPLY)
    {
        exportPLY(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatSTL)
    {
        exportSTL(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatOFF)
    {
        exportOFF(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    LOG_STATUS_IMPORTANT("Exporting Volume Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::exportVolumeGridToMesh(const std::string &prefix,
                                    const bool &formatOBJ, const bool &formatPLY,
                                    const bool &formatOFF, const bool &formatSTL) const
{
    if (!(formatOBJ || formatPLY || formatOFF || formatSTL))
    {
        LOG_WARNING("Exporto mesh option must be enabled to export this mesh. "
                    "User one of the following: "
                    "[--export-obj, --export-ply, --export-off, --export-stl]");
        return;
    }

    TIMER_SET;
    LOG_TITLE("Constructing Volume Grid Mesh");

    // The generated mesh from the volume
    std::unique_ptr< VolumeMesh > volumeGridMesh = std::make_unique< VolumeMesh >();

    // Delta value
    const Vector3f delta(1, 1, 1);

    LOOP_STARTS("Iterating over the grid")
    for (int64_t i = 0; i <  _grid->getWidth(); ++i)
    {
        LOOP_PROGRESS(i, _grid->getWidth());

        for (int64_t j = 0; j < _grid->getHeight(); ++j)
        {
            for (int64_t k = 0; k < _grid->getDepth(); ++k)
            {
                Vector3f coordinate(i, j, k);
                Vector3f pMin = _baseResolution * (coordinate - 0.5f * delta) + _pMin;
                Vector3f pMax = _baseResolution * (coordinate + 0.5f * delta) + _pMin;

                // A mesh representing the bounding box of the cube
                VolumeMesh* voxelCube = VolumeMesh::constructVoxelCube(pMin, pMax);

                // Append it to the volume mesh
                volumeGridMesh->append(voxelCube);

                // Free the voxel cube
                voxelCube->~VolumeMesh();
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOG_STATUS_IMPORTANT("Volume Grid Mesh Construction Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    /// Adjust the dimensions of the resulting mesh
    // Compute the dimensions of the resulting volume mesh
    Vector3f inputMeshbounds = _pMax - _pMin;
    Vector3f inputMeshCenter = _pMin + inputMeshbounds * 0.5;
    Vector3f pMin, pMax, volumeMeshBounds;
    volumeGridMesh->computeBoundingBox(pMin, pMax);
    volumeMeshBounds = pMax - pMin;

    // Compute the scale factor
    Vector3f scaleFactor = inputMeshbounds / volumeMeshBounds;

    // Transform the resulting volume mesh to align with the input mesh
    volumeGridMesh->centerAtOrigin();
    volumeGridMesh->scale(scaleFactor);
    volumeGridMesh->translate(inputMeshCenter);

    LOG_TITLE("Exporting Volume Grid Mesh");
    TIMER_RESET;
    const std::string outputPrefix = prefix + VOLUME_GRID_MESH_SUFFIX;
    if (formatOBJ)
    {
        exportOBJ(outputPrefix,
                  volumeGridMesh->vertices.data(), volumeGridMesh->vertices.size(),
                  volumeGridMesh->triangles.data(), volumeGridMesh->triangles.size());
    }

    if (formatPLY)
    {
        exportPLY(outputPrefix,
                  volumeGridMesh->vertices.data(), volumeGridMesh->vertices.size(),
                  volumeGridMesh->triangles.data(), volumeGridMesh->triangles.size());
    }

    if (formatSTL)
    {
        exportSTL(outputPrefix,
                  volumeGridMesh->vertices.data(), volumeGridMesh->vertices.size(),
                  volumeGridMesh->triangles.data(), volumeGridMesh->triangles.size());
    }

    if (formatOFF)
    {
        exportOFF(outputPrefix,
                  volumeGridMesh->vertices.data(), volumeGridMesh->vertices.size(),
                  volumeGridMesh->triangles.data(), volumeGridMesh->triangles.size());
    }

    LOG_STATUS_IMPORTANT("Exporting Volume Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::exportBoundingBoxMesh(const std::string &prefix,
                                   const bool &formatOBJ, const bool &formatPLY,
                                   const bool &formatOFF, const bool &formatSTL) const
{
    if (!(formatOBJ || formatPLY || formatOFF || formatSTL))
    {
        LOG_WARNING("Exporto mesh option must be enabled to export this mesh. "
                    "User one of the following: "
                    "[--export-obj, --export-ply, --export-off, --export-stl]");
        return;
    }

    TIMER_SET;
    LOG_TITLE("Constructing Volume Bounding Box Mesh");

    // The generated mesh from the volume
    std::unique_ptr< VolumeMesh > volumeMesh = std::make_unique< VolumeMesh >();

    // Delta value
    const Vector3f delta(1, 1, 1);

    // A mesh representing the bounding box of the cube
    VolumeMesh* voxelCube = VolumeMesh::constructVoxelCube(_pMin, _pMax);

    // Append it to the volume mesh
    volumeMesh->append(voxelCube);

    // Free the voxel cube
    delete voxelCube;

    LOG_STATUS_IMPORTANT("Volume Mesh Construction Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    /// Adjust the dimensions of the resulting mesh
    // Compute the dimensions of the resulting volume mesh
    Vector3f inputMeshbounds = _pMax - _pMin;
    Vector3f inputMeshCenter = _pMin + inputMeshbounds * 0.5;
    Vector3f pMin, pMax, volumeMeshBounds;
    volumeMesh->computeBoundingBox(pMin, pMax);
    volumeMeshBounds = pMax - pMin;

    // Compute the scale factor
    Vector3f scaleFactor = inputMeshbounds / volumeMeshBounds;

    // Transform the resulting volume mesh to align with the input mesh
    volumeMesh->centerAtOrigin();
    volumeMesh->scale(scaleFactor);
    volumeMesh->translate(inputMeshCenter);

    LOG_TITLE("Exporting Volume Bounding Box Mesh");
    TIMER_RESET;
    const std::string outputPrefix = prefix + VOLUME_BOUNDING_BOX_MESH_SUFFIX;
    if (formatOBJ)
    {
        exportOBJ(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatPLY)
    {
        exportPLY(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatSTL)
    {
        exportSTL(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatOFF)
    {
        exportOFF(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    LOG_STATUS_IMPORTANT("Exporting Volume Bounding Box Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

uint8_t Volume::getByte(const size_t index) const
{
    return _grid->getByte(index);
}

uint8_t Volume::getConfirmedValue(const int64_t &x,
                                  const int64_t &y,
                                  const int64_t &z) const
{
    if(x > getWidth() - 1 || x < 0 || y > getHeight() - 1 || y < 0 || z > getDepth() -1  || z < 0)
        return 0;

    if (_grid->isFilled(x, y, z))
        return 255;
    else
        return 0;
}

uint8_t Volume::getValueUI8(const int64_t &x,
                            const int64_t &y,
                            const int64_t &z) const
{
    return _grid->getValueUI8(x, y, z);
}

uint16_t Volume::getValueUI16(const int64_t &x,
                              const int64_t &y,
                              const int64_t &z) const
{
    return _grid->getValueUI16(x, y, z);
}

uint32_t Volume::getValueUI32(const int64_t &x,
                              const int64_t &y,
                              const int64_t &z) const
{
    return _grid->getValueUI32(x, y, z);
}

uint64_t Volume::getValueUI64(const int64_t &x,
                              const int64_t &y,
                              const int64_t &z) const
{
    return _grid->getValueUI64(x, y, z);
}

float Volume::getValueF32(const int64_t &x,
                          const int64_t &y,
                          const int64_t &z) const
{
    return _grid->getValueF32(x, y, z);
}

double Volume::getValueF64(const int64_t &x,
                           const int64_t &y,
                           const int64_t &z) const
{
    return _grid->getValueF64(x, y, z);

}

uint8_t Volume::getByte(const int64_t &x,
                        const int64_t &y,
                        const int64_t &z) const
{
    bool outlier;
    uint64_t index = mapToIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return _grid->getByte(index);
}

void Volume::fillVoxel(const int64_t &x,
                       const int64_t &y,
                       const int64_t &z)
{

    if (x - 1 < 0 || x > _grid->getWidth())
    {
        return;
    }
    if (y - 1 < 0 || y > _grid->getHeight())
    {
        return;
    }
    if (z - 1 < 0 || z > _grid->getDepth())
    {
        return;
    }

    _grid->fillVoxel(x - 1, y, z);
    _grid->fillVoxel(x - 1, y - 1, z);
    _grid->fillVoxel(x - 1, y + 1, z);
    _grid->fillVoxel(x - 1, y, z - 1);
    _grid->fillVoxel(x - 1, y - 1, z - 1);
    _grid->fillVoxel(x - 1, y + 1, z - 1);
    _grid->fillVoxel(x - 1, y, z + 1);
    _grid->fillVoxel(x - 1, y - 1, z + 1);
    _grid->fillVoxel(x - 1, y + 1, z + 1);

    // x = 0
    _grid->fillVoxel(x, y - 1, z);
    _grid->fillVoxel(x, y + 1, z);
    _grid->fillVoxel(x, y, z - 1);
    _grid->fillVoxel(x, y - 1, z - 1);
    _grid->fillVoxel(x, y + 1, z - 1);
    _grid->fillVoxel(x, y, z + 1);
    _grid->fillVoxel(x, y - 1, z + 1);
    _grid->fillVoxel(x, y + 1, z + 1);

    // x = +1
    _grid->fillVoxel(x + 1, y, z);
    _grid->fillVoxel(x + 1, y - 1, z);
    _grid->fillVoxel(x + 1, y + 1, z);
    _grid->fillVoxel(x + 1, y, z - 1);
    _grid->fillVoxel(x + 1, y - 1, z - 1);
    _grid->fillVoxel(x + 1, y + 1, z - 1);
    _grid->fillVoxel(x + 1, y, z + 1);
    _grid->fillVoxel(x + 1, y - 1, z + 1);
    _grid->fillVoxel(x + 1, y + 1, z + 1);

}

void Volume::addByte(const size_t &index, const uint8_t byte)
{
    _grid->addByte(index, byte);
}

void Volume::clear(void)
{
    _grid->clear();
}

void Volume::fill(const int64_t &x,
                  const int64_t &y,
                  const int64_t &z)
{
    _grid->fillVoxel(x, y, z);
}

void  Volume::fill(const u_int64_t& index)
{
    _grid->fillVoxel(index);
}

void Volume::clear(const int64_t &x,
                   const int64_t &y,
                   const int64_t &z)
{
    _grid->clearVoxel(x, y, z);
}

void  Volume::clear(const u_int64_t& index)
{
    _grid->clearVoxel(index);
}

size_t Volume::computeNumberNonZeroVoxels(void) const
{
    return _grid->computeNumberNonZeroVoxels();
}

std::string Volume::getFormatString() const
{
    return _grid->getTypeString(_gridType);
}

float Volume::computeVolume() const
{
    // Compute the number of non zero voxels
    const uint64_t numberNonZeroVoxels = _grid->computeNumberNonZeroVoxels();

    // Get the voxel volume in units3
    const float voxelVolume =
            _voxelSize * _voxelSize * _voxelSize;

    // Return the result
    return voxelVolume * numberNonZeroVoxels;
}

void Volume::printStats(const std::string &reference, const std::string *prefix) const
{
    LOG_TITLE("Volume Statistics");

    LOG_STATUS("Collecting Stats.");
    Vector3f bounds = _pMax - _pMin;
    float volumeSize = computeVolume();

    // Write the statistics to a file
    if (prefix != nullptr)
    {
        // Create the file
        std::string fileName = *prefix + "-" + reference + VOLUME_INFO_EXTENSION;
        LOG_STATUS("Writing Info. [ %s ] \n", fileName.c_str());

        FILE* info = fopen(fileName.c_str(), "w");
        fprintf(info, "Stats. [ %s ] \n", reference.c_str());

        if (bounds.x() > 0.f || bounds.y() > 0.f || bounds.z() > 0.f)
        {
            fprintf(info, "\t* Bounding Box:         | [%f, %f, %f] \n",
                     F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
            fprintf(info, "\t* pMin:                 | [%f, %f, %f] \n",
                     F2D(_pMin.x()), F2D(_pMin.y()), F2D(_pMin.z()));
            fprintf(info, "\t* pMax:                 | [%f, %f, %f] \n",
                     F2D(_pMax.x()), F2D(_pMax.y()), F2D(_pMax.z()));
        }

        fprintf(info, "\t* Resolution            | [%d] x [%d] x [%d] \n",
                 I2I32(getWidth()), I2I32(getHeight()), I2I32(getDepth()));
        fprintf(info, "\t* Number of Voxels      | %zu\n",
                getNumberVoxels());
        fprintf(info, "\t* Volume Format         | %s \n",
                 getFormatString().c_str());
        fprintf(info, "\t* Size in Memory        | %sBytes \n",
                 FORMAT(getNumberBytes()));
        fprintf(info, "\t* Volume                | %f³ \n",
                 F2D(volumeSize));

        // Close the file
        fclose(info);
    }

    LOG_STATUS("Volume [ %s ]", reference.c_str());

    if (bounds.x() > 0.f || bounds.y() > 0.f || bounds.z() > 0.f)
    {
        LOG_INFO("\t* Bounding Box:         | [%f, %f, %f]",
                 F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
        LOG_INFO("\t* pMin:                 | [%f, %f, %f]",
                 F2D(_pMin.x()), F2D(_pMin.y()), F2D(_pMin.z()));
        LOG_INFO("\t* pMax:                 | [%f, %f, %f]",
                 F2D(_pMax.x()), F2D(_pMax.y()), F2D(_pMax.z()));
    }

    LOG_INFO("\t* Resolution            | [%d] x [%d] x [%d]",
             I2I32(getWidth()), I2I32(getHeight()), I2I32(getDepth()));
    LOG_INFO("\t* Number of Voxels      | %" PRIu64 "",
             getNumberVoxels());
    LOG_INFO("\t* Volume Format         | %s",
             getFormatString().c_str());
    LOG_INFO("\t* Size in Memory        | %sBytes",
             FORMAT(getNumberBytes()));
    LOG_INFO("\t* Volume                | %f³",
             F2D(volumeSize));
}

void Volume::addVolumePass(const Volume* volume)
{
    // Start the timer
    TIMER_SET;

    // If the volume dimensions are not similar to this one
    // then, print an ERROR ...
    if (!(getWidth() == volume->getWidth() &&
          getHeight() == volume->getHeight() &&
          getDepth() == volume->getDepth()))
    {
        LOG_ERROR("The dimensions of the two volumes don't match");
        return;
    }

    // Loop over the volume elements byte-by-byte and add them to the corresponding one in this voume
    for (size_t i = 0; i < volume->getNumberBytes(); ++i)
    {
        addByte(i, volume->getByte(i));
    }

    _addingVolumePassTime = GET_TIME_SECONDS;
}

void Volume::addVolume(const std::string &volumePrefix)
{
//    Volume* volume = new Volume(volumePrefix, VOLUME_TYPE::UI8);

//    // If the volume dimensions are not similar to this one
//    // then, print an ERROR ...
//    if (!(getWidth() == volume->getWidth() &&
//          getHeight() == volume->getHeight() &&
//          getDepth() == volume->getDepth()))
//    {
//        LOG_ERROR("The dimensions of the two volumes don't match");
//        return;
//    }

//    // Loop over the volume elements byte-by-byte and add
//    // them to the corresponding one in this voume
//    for (size_t i = 0; i < volume->getNumberBytes(); ++i)
//    {
//        addByte(i, volume->getByte(i));
//    }

//    // Release the gird to save some memory
//    delete volume;
}

size_t Volume::getNumberBytes(void) const
{
    return _grid->getNumberBytes();
}

int32_t Volume::getLargestDimension(const Vector3f& dimensions)
{
    uint32_t index = 0;
    float value = dimensions[index];


    for (int32_t i = 1; i < DIMENSIONS; ++i)
    {
        if (value < dimensions[i])
        {
            index = i;
            value = dimensions[i];
        }
    }

    return index;
}

Volume::~Volume()
{
    delete _grid;
}


Volume* Volume::constructIsoValueVolume(const Volume* volume,
                                        const size_t& isoValue)
{
    // Create the iso-volume
    Volume* isoVolume = new Volume(volume->getWidth(), volume->getHeight(),volume->getDepth(),
                                   volume->getPMin(), volume->getPMax(),
                                   VOLUME_TYPE::BIT);

    LOG_STATUS("Constructing Iso Volume from a Single Value [%zu]", isoValue);
    for (size_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                if (volume->getValueUI64(x, y, z) == isoValue)
                {
                    isoVolume->fill(x, y, z);
                }
                else
                {
                    isoVolume->clear(x, y, z);
                }
            }
        }
    }
    LOOP_DONE;

    // Return the resulting iso-volume
    return isoVolume;
}

Volume* Volume::constructVolumeWithMinimumIsoValue(const Volume* volume,
                                                   const size_t& minIsoValue)
{
    // Create the iso-volume
    Volume* isoVolume = new Volume(volume->getWidth(), volume->getHeight(),volume->getDepth(),
                                   volume->getPMin(), volume->getPMax(),
                                   VOLUME_TYPE::BIT);

    LOG_STATUS("Constructing Iso Volume with Minimum Value [%zu]", minIsoValue);
    for (size_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                if (volume->getValueUI64(x, y, z) >= minIsoValue)
                {
                    isoVolume->fill(x, y, z);
                }
                else
                {
                    isoVolume->clear(x, y, z);
                }
            }
        }
    }
    LOOP_DONE;

    // Return the resulting iso-volume
    return isoVolume;
}

Volume* Volume::constructVolumeWithMaximumIsoValue(const Volume* volume,
                                                   const size_t& maxIsoValue)
{
    // Create the iso-volume
    Volume* isoVolume = new Volume(volume->getWidth(), volume->getHeight(),volume->getDepth(),
                                   volume->getPMin(), volume->getPMax(),
                                   VOLUME_TYPE::BIT);

    LOG_STATUS("Constructing Iso Volume with Minimum Value [%zu]", maxIsoValue);
    for (size_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                if (volume->getValueUI64(x, y, z) <= maxIsoValue)
                {
                    isoVolume->fill(x, y, z);
                }
                else
                {
                    isoVolume->clear(x, y, z);
                }
            }
        }
    }
    LOOP_DONE;

    // Return the resulting iso-volume
    return isoVolume;
}

Volume* Volume::constructVolumeWithIsoRange(const Volume* volume,
                                            const size_t& minIsoValue,
                                            const size_t& maxIsoValue)
{
    // Create the iso-volume
    Volume* isoVolume = new Volume(volume->getWidth(), volume->getHeight(),volume->getDepth(),
                                   volume->getPMin(), volume->getPMax(),
                                   VOLUME_TYPE::BIT);

    LOG_STATUS("Constructing Iso Volume with Range [%zu - %zu]", minIsoValue, maxIsoValue);
    for (size_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                if (volume->getValueUI64(x, y, z) <= maxIsoValue &&
                        volume->getValueUI64(x, y, z) >= minIsoValue)
                {
                    isoVolume->fill(x, y, z);
                }
                else
                {
                    isoVolume->clear(x, y, z);
                }
            }
        }
    }
    LOOP_DONE;

    // Return the resulting iso-volume
    return isoVolume;
}

Volume* Volume::constructNonZeroVolume(const Volume* volume)
{
    // Create the iso-volume
    Volume* isoVolume = new Volume(volume->getWidth(), volume->getHeight(),volume->getDepth(),
                                   volume->getPMin(), volume->getPMax(),
                                   VOLUME_TYPE::BIT);

    LOG_STATUS("Constructing Full Range Iso Volume");
    for (size_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                if (volume->getValueUI64(x, y, z) > 0)
                {
                    isoVolume->fill(x, y, z);
                }
                else
                {
                    isoVolume->clear(x, y, z);
                }
            }
        }
    }
    LOOP_DONE;

    // Return the resulting iso-volume
    return isoVolume;
}

Volume* Volume::constructIsoValuesVolume(const Volume* volume,
                                         const std::vector<size_t> &isoValues)
{
    // Create the iso-volume
    Volume* isoVolume = new Volume(volume->getWidth(), volume->getHeight(),volume->getDepth(),
                                   volume->getPMin(), volume->getPMax(),
                                   VOLUME_TYPE::BIT);

    LOG_STATUS("Constructing Iso Volume");
    for (int64_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                for (size_t iv = 0; iv < isoValues.size(); ++iv)
                {
                    const auto &voxelValue = volume->getValueUI64(x, y, z);
                    if (std::find(isoValues.begin(), isoValues.end(), voxelValue) != isoValues.end())
                    {
                        isoVolume->fill(x, y, z);
                    }
                    else
                    {
                        isoVolume->clear(x, y, z);
                    }
                }
            }
        }
    }
    LOOP_DONE;

    return isoVolume;
}

std::vector<size_t> Volume::createHistogram(const Volume* volume,
                                            const VOLUME_TYPE& type)
{
    uint64_t histogramWidth;
    switch (type)
    {
    case VOLUME_TYPE::UI8:
    {
        histogramWidth = std::numeric_limits<uint8_t>::max();
    } break;

    case VOLUME_TYPE::UI16:
    {
        histogramWidth = std::numeric_limits<uint16_t>::max();
    } break;

    case VOLUME_TYPE::UI32:
    {
        histogramWidth = std::numeric_limits<uint32_t>::max();
    } break;

    case VOLUME_TYPE::UI64:
    {
        histogramWidth = std::numeric_limits<uint64_t>::max();
    } break;

    case VOLUME_TYPE::BIT:
    {
        LOG_ERROR("Histograms CANNOT be computed to a bit volume!");
    }

    case VOLUME_TYPE::F32:
    case VOLUME_TYPE::F64:
    default:
        LOG_ERROR("Histograms CANNOT be computed to float volumes!");
        break;
    }

    // Allocation
    std::vector< size_t > histogram;
    histogram.resize(histogramWidth);

    // Initialization
    for (size_t i = 0; i < histogramWidth; ++i)
        histogram[i] = 0;

    for (size_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (size_t y = 0; y < volume->getHeight(); ++y)
        {
            for (size_t z = 0; z < volume->getDepth(); ++z)
            {
                histogram[volume->getValueUI64(x, y, z)] += 1;
            }
        }
    }

    // Return the histogram array
    return histogram;
}

Volume* Volume::constructFromTiffMask(
        const std::string &maskDirectory,
        const int64_t &maskWidth, const int64_t &maskHeight,
        const Ultraliser::VOLUME_TYPE& gridType)
{
#ifdef ULTRALISER_USE_TIFF
    // Set the timer
    TIMER_SET;

    LOG_TITLE("Loading .Tiff");
    LOG_SUCCESS("Mask Directory [ %s ]", maskDirectory.c_str());

    // Get a list of all the stacks of the masks
    std::vector< std::string > maskFiles;
    Ultraliser::Directory::list(maskDirectory, maskFiles, ".tif");

    // Sort the mask files
    std::sort(maskFiles.begin(), maskFiles.end());

    // Adding a little delta
    int64_t numZeroPaddingVoxels = 16;
    Volume* maskVolume = new Volume(
                maskWidth + numZeroPaddingVoxels,
                maskHeight + numZeroPaddingVoxels,
                I2I64(maskFiles.size()) + numZeroPaddingVoxels,
                Vector3f(),
                Vector3f(),
                gridType);
    LOG_SUCCESS("Mask Dimensions [%d x %d x %d]",
                maskVolume->getWidth(), maskVolume->getHeight(), maskVolume->getDepth());

    PROGRESS_SET;
    LOG_STATUS("Constructing Volume");
    OMP_PARALLEL_FOR
    for(int64_t i = 0; i < I2I64(maskFiles.size()); ++i)
    {
        // Read the image
        std::string imagePath = maskDirectory + "/" + maskFiles[I2UI64(i)];
        std::unique_ptr< TiffImage > image(new TiffImage);
        image->setimageFile(imagePath);
        image->readImage();

        // Update the volume
        for(int64_t x = 0; x < maskWidth; ++x)
        {
            for(int64_t y = 0; y < maskHeight; ++y)
            {
                if(image->isPixelFilled(I2I32(x), I2I32(y)))
                {
                    maskVolume->fill(x + numZeroPaddingVoxels / 2,
                                     y + numZeroPaddingVoxels / 2,
                                     i + numZeroPaddingVoxels / 2);
                }
            }
        }

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, maskFiles.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;

    // Statictics
    LOG_STATS(GET_TIME_SECONDS);

    // Return a pointer to the mask volume
    return maskVolume;
#else
    LOG_ERROR("Ultraliser compiled with NO support to read .TIFF stacks!");
    return nullptr;
#endif
}

Volume::SOLID_VOXELIZATION_AXIS Volume::getSolidvoxelizationAxis(const std::string &argumentString)
{
    if (argumentString == "x")
    {
        return SOLID_VOXELIZATION_AXIS::X;
    }
    else if (argumentString == "y")
    {
        return SOLID_VOXELIZATION_AXIS::Y;
    }
    else if (argumentString == "z")
    {
        return SOLID_VOXELIZATION_AXIS::Z;
    }
    else if (argumentString == "xyz")
    {
        return SOLID_VOXELIZATION_AXIS::XYZ;
    }
    else
    {
        // Error, therefore terminate
        LOG_ERROR("The option [ %s ] is not valid for --solid-voxelization-axis! "
                  "Please use one of the following [x, y, z, xyz].", argumentString.c_str());

        // For the sake of compilation only.
        return SOLID_VOXELIZATION_AXIS::XYZ;
    }
}

bool Volume::isBorderVoxel(const int64_t& x, const int64_t& y,const int64_t& z) const
{
    if(!isFilledWithoutBoundCheck(x, y, z))
        return false;

    if(!isFilledWithoutBoundCheck(x + 1, y, z))
        return true;

    if(!isFilledWithoutBoundCheck(x - 1, y, z))
        return true;

    if(!isFilledWithoutBoundCheck(x, y + 1, z))
        return true;

    if(!isFilledWithoutBoundCheck(x, y - 1, z))
        return true;

    if(!isFilledWithoutBoundCheck(x, y, z + 1))
        return true;

    if(!isFilledWithoutBoundCheck(x, y, z - 1))
        return true;

    return false;
}

std::vector< CandidateVoxels > Volume::searchForCandidateVoxels() const
{
    // This list will collect the border voxels per slice (along the width)
    std::vector< CandidateVoxels > perSliceBorderVoxels;
    perSliceBorderVoxels.reserve(getWidth());

    // Collect the border voxels per slice in parallel
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < getWidth(); ++i)
    {
        auto& borderVoxel = perSliceBorderVoxels[i];

        for (size_t j = 0; j < getHeight(); ++j)
        {
            for (size_t k = 0; k < getDepth(); ++k)
            {
                if (isBorderVoxel(i, j, k))
                {
                    CandidateVoxel* voxel = new CandidateVoxel();
                    voxel->x = i;
                    voxel->y = j;
                    voxel->z = k;
                    voxel->deletable = false;
                    borderVoxel.push_back(voxel);
                }
            }
        }
    }

    return perSliceBorderVoxels;
}


#include <algorithm>
size_t Volume::deleteCandidateVoxels(std::unique_ptr< Thinning6Iterations > &thinning)
{
    size_t numberDeletedVoxels = 0;

    // This list will collect the border voxels per slice (along the width)
    std::vector< CandidateVoxels > perSliceBorderVoxels;
    perSliceBorderVoxels.resize(getWidth());

    // Collect the border voxels per slice in parallel

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < getWidth(); ++i)
    {
        for (size_t j = 0; j < getHeight(); ++j)
        {
            for (size_t k = 0; k < getDepth(); ++k)
            {
                if (isBorderVoxel(i, j, k))
                {
                    CandidateVoxel* voxel = new CandidateVoxel();
                    voxel->x = i;
                    voxel->y = j;
                    voxel->z = k;
                    voxel->deletable = false;

                    perSliceBorderVoxels[i].push_back(voxel);
                }

            }
        }
    }

    for (size_t direction = 0; direction < 6; direction++)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
        {
             size_t j, k;
            for (j = 0; j < perSliceBorderVoxels[i].size(); ++j)
            {
                // A block of the volume that is scanned every iteration
                int8_t volumeBlock[26];

                for (k = 0; k < 26; k++)
                {
                    size_t idx, idy, idz;

                    idx = perSliceBorderVoxels[i][j]->x + VDX[k];
                    idy = perSliceBorderVoxels[i][j]->y + VDY[k];
                    idz = perSliceBorderVoxels[i][j]->z + VDZ[k];

                    volumeBlock[k] = isFilledWithoutBoundCheck(idx, idy, idz) ? 1 : 0;
                }

                if (thinning->matches(direction, volumeBlock))
                {
                    perSliceBorderVoxels[i][j]->deletable = true;
                }
            }
        }

        for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
        {
            for (size_t j = 0; j < perSliceBorderVoxels[i].size(); ++j)
            {
                if (perSliceBorderVoxels[i][j]->deletable)
                {
                    numberDeletedVoxels++;
                    clear(perSliceBorderVoxels[i][j]->x,
                          perSliceBorderVoxels[i][j]->y,
                          perSliceBorderVoxels[i][j]->z);

                    perSliceBorderVoxels[i][j]->deletable = false;
                }
            }
        }
    }

    // Clear the border voxels (list of lists)
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
    {
        perSliceBorderVoxels[i].clear();
        perSliceBorderVoxels[i].shrink_to_fit();
    }

    perSliceBorderVoxels.clear();
    perSliceBorderVoxels.shrink_to_fit();


    return numberDeletedVoxels;
}

CandidateVoxels Volume::searchForCandidateVoxelsOne() const
{
    // This list will collect the border voxels per slice (along the width)
    std::vector< CandidateVoxels > perSliceBorderVoxels;
    perSliceBorderVoxels.resize(getWidth());

    // Collect the border voxels per slice in parallel
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < getWidth(); ++i)
    {
        for (size_t j = 0; j < getHeight(); ++j)
        {
            for (size_t k = 0; k < getDepth(); ++k)
            {
                if (isBorderVoxel(i, j, k))
                {
                    CandidateVoxel* voxel = new CandidateVoxel();
                    voxel->x = i;
                    voxel->y = j;
                    voxel->z = k;
                    voxel->deletable = false;

                    perSliceBorderVoxels[i].push_back(voxel);
                }

            }
        }
    }

    size_t allSize = 0;
    for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
    {
        allSize += perSliceBorderVoxels[i].size();
    }

    CandidateVoxels candiateVoxels;
    candiateVoxels.reserve(allSize);
    for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
    {
        candiateVoxels.insert(candiateVoxels.end(),
                              perSliceBorderVoxels[i].begin(),
                              perSliceBorderVoxels[i].end());

        perSliceBorderVoxels[i].clear();
        perSliceBorderVoxels[i].shrink_to_fit();
    }

    return candiateVoxels;
}


std::vector< std::vector< Vec3ui_64 > > Volume::searchForBorderVoxels() const
{
    // This list will collect the border voxels per slice (along the width)
    std::vector< std::vector< Vec3ui_64 > > perSliceBorderVoxels;
    perSliceBorderVoxels.resize(getWidth());

    // Collect the border voxels per slice in parallel
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < getWidth(); ++i)
    {
        for (size_t j = 0; j < getHeight(); ++j)
        {
            for (size_t k = 0; k < getDepth(); ++k)
            {
                if (isBorderVoxel(i, j, k))
                {
                    perSliceBorderVoxels[i].push_back(Vec3ui_64(i, j, k));
                }
            }
        }
    }

    return perSliceBorderVoxels;
}


void Volume::confirmDeletableVoxels(CandidateVoxels& candidateVoxels,
                                    std::unique_ptr< Thinning6Iterations > &thinning,
                                    int direction) const
{
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < candidateVoxels.size(); ++i)
    {
        // A block of the volume that is scanned every iteration
        int8_t volumeBlock[26];

        for (size_t k = 0; k < 26; k++)
        {
            size_t idx, idy, idz;

            idx = candidateVoxels[i]->x + VDX[k];
            idy = candidateVoxels[i]->y + VDY[k];
            idz = candidateVoxels[i]->z + VDZ[k];

            volumeBlock[k] = isFilledWithoutBoundCheck(idx, idy, idz) ? 1 : 0;
        }

        if (thinning->matches(direction, volumeBlock))
        {
            candidateVoxels[i]->deletable = true;
        }
    }
}

std::vector< Vec3ui_64 >
Volume::searchForDeletableVoxels(std::vector< std::vector< Vec3ui_64 > > &perSliceBorderVoxels,
                                 std::unique_ptr< Thinning6Iterations > &thinning,
                                 int direction) const
{
    std::vector< Vec3ui_64 > voxelsToBeDeleted;


    for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
    {
        for (size_t j = 0; j < perSliceBorderVoxels[i].size(); ++j)
        {
            const auto& borderVoxel = perSliceBorderVoxels[i][j];

            // A block of the volume that is scanned every iteration
            int8_t volumeBlock[26];

            //
            for (size_t k = 0; k < 26; k++)
            {
                const auto idx = borderVoxel[0] + VDX[k];
                const auto idy = borderVoxel[1] + VDY[k];
                const auto idz = borderVoxel[2] + VDZ[k];

                volumeBlock[k] = isFilledWithoutBoundCheck(idx, idy, idz) ? 1 : 0;
            }

            if (thinning->matches(direction, volumeBlock))
            {
                voxelsToBeDeleted.push_back(perSliceBorderVoxels[i][j]);
            }
        }
    }
    return voxelsToBeDeleted;
}

















struct GraphEdge
{
    GraphEdge() {}
    GraphEdge(size_t index, size_t index0, size_t index1)
    {
        i = index;
        p0Index = index0;
        p1Index = index1;
    }

    size_t i;
    size_t p0Index;
    size_t p1Index;
};

struct GraphNode
{
    GraphNode() { }
    GraphNode(Vector3f point, Vector3f voxel, size_t index)
    {
        this->pNode = point;
        this->pVoxel = voxel;
        this->index = index;
    }

    Vector3f pNode;
    Vector3f pVoxel;

    size_t index = 0;

    float radius;

    // If the node has been visited before
    bool visited = false;

    // If the node represents a branching point
    bool branching = false;

    // If the node is a terminal one, i.e. only has a single edge
    bool terminal = false;

    // If the node is located inside the soma
    bool insideSoma = false;

    bool isSoma = false;

    // If the node is connected to soma (node of emanating branch from the soma)
    bool connectedToSoma = false;

    std::vector< GraphNode* > edgeNodes;
};


struct Branch
{
    // A list of points
    std::vector< GraphNode* > nodes;

    size_t index;

    bool valid = true;

    Branch* parent = nullptr;
    std::vector< Branch* > children;
};


typedef std::vector< Branch* > Branches;




Branch* buildBranchTillNextBranchingPoint(GraphNode* branchingNode, GraphNode* edgeNode)
{
    Branch* branch = new Branch();

    branch->nodes.push_back(branchingNode);
    branch->nodes.push_back(edgeNode);

    if (edgeNode->terminal)
        return branch;

    if (edgeNode->branching)
        return branch;

    GraphNode *aux0 = branchingNode;
    GraphNode *aux1 = edgeNode;

    while (aux1->edgeNodes.size() == 2)
    {
        auto n0 = aux1->edgeNodes[0];
        auto n1 = aux1->edgeNodes[1];

        if (n0->index == aux0->index)
        {
            aux0 = aux1;
            aux1 = n1;
        }
        else
        {
            aux0 = aux1;
            aux1 = n0;
        }

        branch->nodes.push_back(aux1);
    }

    return branch;
}

Branch* buildBranch(GraphNode* firstNode, GraphNode* edgeNode)
{
    Branch* branch = new Branch();

    firstNode->visited = true;
    edgeNode->visited = true;

    branch->nodes.push_back(firstNode);
    branch->nodes.push_back(edgeNode);

    // If the node along the edge is a terminal, then append this branch to the list and return
    if (edgeNode->edgeNodes.size() == 1)
    {
        return branch;
    }

    // If the node along the edge is a branching node, then append this branch to the list and
    // consider the next branch
    else if (edgeNode->edgeNodes.size() > 2)
    {
        return branch;
    }

    // Otherwise, construct the branches from connected segments
    else
    {
        // The previous node is the first node
        GraphNode *previousNode = firstNode;

        // The current node is the edge node
        GraphNode *currentNode = edgeNode;

        // Ensure that the current node has only two connected edges (or nodes)
        while (true)
        {
            // Get a reference to the connecting nodes to the current node
            auto edgeNode0 = currentNode->edgeNodes[0];
            auto edgeNode1 = currentNode->edgeNodes[1];

            // Ignore the previous node
            if (edgeNode0->index == previousNode->index)
            {
                previousNode = currentNode;
                currentNode = edgeNode1;
            }
            else
            {
                previousNode = currentNode;
                currentNode = edgeNode0;
            }

            branch->nodes.push_back(currentNode);
            currentNode->visited = true;

            if (!(currentNode->edgeNodes.size() == 2))
                break;
        }
    }

    return branch;

}
Branches buildBranchesFromNodes(std::vector< GraphNode* > nodes)
{
    Branches branches;
    size_t branchIndex = 0;

    // Construct the heirarichy to the terminal
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        auto& node = nodes[i];

        // The node must be branching
        if (node->branching)
        {
            // Construct the branch
            for (size_t j = 0; j < node->edgeNodes.size(); ++j)
            {
                auto& edgeNode = node->edgeNodes[j];
                {
                    // If the edgeNode is visited before, then this branch has been reconstructed
                    if (edgeNode->visited)
                        continue;

                    auto branch = buildBranch(node, edgeNode);
                    branch->index = branchIndex;
                    branchIndex++;

                    branches.push_back(branch);
                }
            }
        }
    }

    return branches;
}

Volume* Volume::getBrick(const size_t& x1, const size_t& x2,
                         const size_t& y1, const size_t& y2,
                         const size_t& z1, const size_t& z2)
{
    // Ensure that the given grid coordinates are correct
    const auto minX = x1 < x2 ? x1 : x2;
    const auto maxX = x1 > x2 ? x1 : x2;

    const auto minY = y1 < y2 ? y1 : y2;
    const auto maxY = y1 > y2 ? y1 : y2;

    const auto minZ = z1 < z2 ? z1 : z2;
    const auto maxZ = z1 > z2 ? z1 : z2;

    // Determine the correct dimensions
    const int64_t width = maxX - minX;
    const int64_t height = maxY - minY;
    const int64_t depth = maxZ - minZ;

    // Adjust the pMin and pMax relatively

    // Create the volume
    Volume* brick = new Volume(width, height, depth);

    for (size_t i = 0; i < width; ++i)
    {
        for (size_t j = 0; j < height; ++j)
        {
            for (size_t k = 0; k < depth; ++k)
            {
                if (isFilledWithoutBoundCheck(i + minX, j + minY, k + minZ))
                {
                    brick->fillVoxel(i, j, k);
                }
                else
                {
                    brick->clear(i, j, k);
                }
            }
        }
    }

    // Return the created volume brick
    return brick;
}

std::vector< Vector3f > Volume::applyThinning(std::vector< Vector3f > & centers, std::vector< float >& radii)
{
    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel = std::make_unique<Thinning6Iterations>();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    // Compute the dimensions of the resulting volume mesh
    const Vector3f inputMeshbounds = _pMax - _pMin;
    const Vector3f inputMeshCenter = _pMin + inputMeshbounds * 0.5f;

    // Vector3f pMinVolume, pMaxVolume, boundsVolume;
    const Vector3f pMaxVolume(getWidth() * 1.f, getHeight() * 1.f, getDepth() * 1.f);
    const Vector3f volumeCenter = (0.5f * pMaxVolume);
    const Vector3f scaleFactor = inputMeshbounds / pMaxVolume;

    std::vector< Vector3f > shellPoints;
    std::vector< std::vector< Vec3ui_64 > > perSliceSurfaceShell = searchForBorderVoxels();
    for (size_t i = 0; i < perSliceSurfaceShell.size(); ++i)
    {
        for (size_t j = 0; j < perSliceSurfaceShell[i].size(); ++j)
        {
            const auto voxel = perSliceSurfaceShell[i][j];
            shellPoints.push_back(Vector3f(voxel.x(), voxel.y(), voxel.z()));
        }
        perSliceSurfaceShell[i].clear();
    }
    perSliceSurfaceShell.clear();

    // Adjust the locations of the shell points
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < shellPoints.size(); ++i)
    {
        // Center at the origin
        shellPoints[i] -= volumeCenter;

        // Scale
        shellPoints[i].x() *= scaleFactor.x();
        shellPoints[i].y() *= scaleFactor.y();
        shellPoints[i].z() *= scaleFactor.z();

        // Translate
        shellPoints[i] += inputMeshCenter;
    }

    TIMER_SET;
    LOG_STATUS("Thinning Volume");
    LOOP_STARTS("Thinning Loop");
    LOOP_PROGRESS(0, 100);
    while(1)
    {
        size_t numberDeletedVoxels = 0;

        // std::vector< std::vector< Vec3ui_64 > > perSliceBorderVoxels = searchForBorderVoxels();

        CandidateVoxels candidateVoxels = searchForCandidateVoxelsOne();
        std::cout << "Hola \n";

        for (size_t direction = 0; direction < 6; direction++)
        {
            // Search for the delerable voxels
            // std::vector< Vec3ui_64 > voxelsToBeDeleted =
            //        searchForDeletableVoxels(
            //            perSliceBorderVoxels, thinningKernel, direction);

            confirmDeletableVoxels(candidateVoxels, thinningKernel, direction);

            // Delete the voxels
            for (size_t i = 0; i < candidateVoxels.size(); ++i)
            {
                if (candidateVoxels[i]->deletable)
                {
                    numberDeletedVoxels++;
                    clear(candidateVoxels[i]->x, candidateVoxels[i]->y, candidateVoxels[i]->z);
                }
            }

            // Clear the container of the deleted voxels
            // voxelsToBeDeleted.clear();
        }

        // Clear all the containers of the border voxels
//        for (size_t i = 0; i < perSliceBorderVoxels.size(); ++i)
//        {
//            perSliceBorderVoxels[i].clear();
//        }
//         perSliceBorderVoxels.clear();

         // Updating the progess bar
        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
        //LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
        //              initialNumberVoxelsToBeDeleted);

        if (numberDeletedVoxels == 0)
            break;

        loopCounter++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);


    std::map< size_t, size_t > volumeToNodeMap;


    // Construct all the nodes
    size_t nodeIndex =0;
    std::vector< GraphNode* > nodes;
    for (size_t i = 0; i < getWidth(); ++i)
    {
        for (size_t j = 0; j < getHeight(); ++j)
        {
            for (size_t k = 0; k < getDepth(); ++k)
            {
                if (isFilled(i, j, k))
                {
                    size_t voxelIndex = mapTo1DIndexWithoutBoundCheck(i, j, k);

                    Vector3f pVoxel(i * 1.f, j * 1.f, k * 1.f);

                    Vector3f pNode(pVoxel);
                    pNode -= volumeCenter;
                    pNode.x() *= scaleFactor.x();
                    pNode.y() *= scaleFactor.y();
                    pNode.z() *= scaleFactor.z();
                    pNode += inputMeshCenter;

                    nodes.push_back(new GraphNode(pNode, pVoxel, voxelIndex));
                    volumeToNodeMap.insert(std::pair< size_t, size_t >(voxelIndex, nodeIndex));
                    nodeIndex++;
                }
            }
        }
    }

    // Calculate the radii of every point
    std::vector< float > nodesRadii;
    nodesRadii.resize(nodes.size());

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        float minimumDistance = 1e32;
        for (size_t j = 0; j < shellPoints.size(); ++j)
        {
            const float distance = (nodes[i]->pNode - shellPoints[j]).abs();

            if (distance < minimumDistance)
                minimumDistance = distance;
        }
        nodesRadii[i] = minimumDistance;
    }


    auto it = std::max_element(std::begin(nodesRadii), std::end(nodesRadii));
    auto largestRadiusIndex = std::distance(std::begin(nodesRadii), it);

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->radius = nodesRadii[i];
    }
    nodesRadii.clear();

    // Rasterize a sphere
    Vector3f somaCenter(nodes[largestRadiusIndex]->pNode.x(),
                        nodes[largestRadiusIndex]->pNode.y(),
                        nodes[largestRadiusIndex]->pNode.z());
    Vector3f somaCenterVoxel(nodes[largestRadiusIndex]->pVoxel.x(),
                             nodes[largestRadiusIndex]->pVoxel.y(),
                             nodes[largestRadiusIndex]->pVoxel.z());
    const auto somaRadius = nodes[largestRadiusIndex]->radius;
    Sample* somaSphere  = new Sample(somaCenter, somaRadius, 0);

    std::cout << "Soma: "
              << somaSphere->getPosition().x() << " "
              << somaSphere->getPosition().y() << " "
              << somaSphere->getPosition().z() << " "
              << somaSphere->getRadius() << "\n";

    // Construct the graph and connect the nodes
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // Check if the node has been visited before
        GraphNode* node = nodes[i];

        // Count the number of the connected edges to the node
        size_t connectedEdges = 0;

        // Search for the neighbours
        for (size_t l = 0; l < 26; l++)
        {
            size_t idx = node->pVoxel.x() + VDX[l];
            size_t idy = node->pVoxel.y() + VDY[l];
            size_t idz = node->pVoxel.z() + VDZ[l];

            if (isFilled(idx, idy, idz))
            {
                connectedEdges++;

                // Find the index of the voxel
                const auto& voxelIndex = mapTo1DIndexWithoutBoundCheck(idx, idy, idz);

                // Find the corresponding index of the node to access the node from the nodes list
                const auto& nodeIndex = volumeToNodeMap.find(voxelIndex)->second;

                // Add the node to the edgeNodes, only to be able to access it later
                node->edgeNodes.push_back(nodes[nodeIndex]);
            }
        }

        if (connectedEdges == 1)
            node->terminal = true;
        else if (connectedEdges > 2)
            node->branching = true;

    }

//    // Verify the nodes that are located within the soma
//    // OMP_PARALLEL_FOR
//    size_t insideSomaSample = 0;
//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        GraphNode* node = nodes[i];

//        if (isPointInSphere(node->pNode, somaCenter, somaRadius))
//        {
//            node->insideSoma = true;
//            insideSomaSample++;
//        }
//    }

    // Build soma node
    GraphNode* somaNode = new GraphNode(somaCenter, somaCenterVoxel, nodeIndex++);
    somaNode->radius = somaRadius;
    somaNode->isSoma = true;
    nodes.push_back(somaNode);

    // Re-index the samples
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->index = i;
    }






//    // OMP_PARALLEL_FOR
//    size_t arbors = 0;
//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        GraphNode* node = nodes[i];

//        if (!node->insideSoma)
//        {
//            // Detect if the node is located on the soma
//            for (size_t j = 0; j < node->edgeNodes.size(); ++j)
//            {
//                if (node->edgeNodes[j]->insideSoma)
//                {
//                    node->connectedToSoma = true;
//                    arbors++;
//                    break;
//                }
//            }

//            // Disconnect the samples that are located in the soma
//            std::vector< GraphNode* > edgeNodes;
//            bool connectToSomaNode = false;
//            for (size_t j = 0; j < node->edgeNodes.size(); ++j)
//            {
//                // If the edge node is not located in the soma, then add it to the edgeNodes list
//                if (!node->edgeNodes[j]->insideSoma)
//                {
//                    edgeNodes.push_back(node->edgeNodes[j]);
//                }
//                else
//                {
//                    // Connect to the soma node flag
//                    connectToSomaNode = true;
//                }
//            }

//            // If connected to soma, then connect to soma
//            if (connectToSomaNode)
//            {
//                edgeNodes.push_back(somaNode);
//            }

//            // Clear the current edgeNodes vector
//            node->edgeNodes.clear();
//            node->edgeNodes.shrink_to_fit();

//            // Connect it to the new vector
//            node->edgeNodes = edgeNodes;
//        }
//    }

    Branches branches = buildBranchesFromNodes(nodes);

    std::cout << "Branches: " << branches.size() << "\n";
    std::cout << "Nodes (Samples): " << nodes.size() << "\n";

    // Filter the branches the are located inside the soma
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < branches.size(); ++i)
    {
        auto& branch = branches[i];

        size_t countSamplesInsideSoma = 0;
        for (size_t j = 0; j < branch->nodes.size(); ++j)
        {
            if (isPointInSphere(branch->nodes[j]->pNode, somaCenter, somaRadius))
            {
                branch->nodes[j]->insideSoma = true;
                countSamplesInsideSoma++;
            }
        }

        // If the count of the sample located inside the soma is zero, then it is a valid branch
        if (countSamplesInsideSoma == 0)
        {
            branch->valid = true;
        }

        // If all the branch nodes are located inside the soma, then it is not valid
        else if (countSamplesInsideSoma == branch->nodes.size())
        {
            branch->valid = false;
        }

        // Otherwise, it is a branch that is connected to the soma
        else
        {
            std::vector< GraphNode* > newNodes;

            // Get the first and last nodes
            auto& firstNode = branch->nodes.front();
            auto& lastNode = branch->nodes.back();

            if (firstNode->insideSoma)
            {
                newNodes.push_back(somaNode);
                for (size_t j = 0; j < branch->nodes.size(); ++j)
                {
                    if (branch->nodes[j]->insideSoma)
                        continue;
                    else
                    {
                        newNodes.push_back(branch->nodes[j]);
                    }
                }
            }
            else if (lastNode->insideSoma)
            {
                for (size_t j = 0; j < branch->nodes.size(); ++j)
                {
                    if (branch->nodes[j]->insideSoma)
                        continue;
                    else
                    {
                        newNodes.push_back(branch->nodes[j]);
                    }
                }
                newNodes.push_back(somaNode);
            }

            branch->nodes.clear();
            branch->nodes.shrink_to_fit();
            branch->nodes = newNodes;
            branch->valid = true;
        }
    }


    // sCenter = somaCenter;
    // sRadius = somaRadius;





//    std::cout << "Writing \n";

    std::fstream xstream;
//    xstream.open("/abdellah2/scratch/thinning/output/projections/nodes.txt", std::ios::out);

//    for (size_t i = 0; i < branches.size(); ++i)
//    {
//        if (branches[i]->nodes.size() == 0)
//        {
//            std::cout << "Empty branch \n";
//            continue;
//        }

//        if (!branches[i]->valid)
//        {
//            std::cout << "Invalid branch \n";
//            continue;
//        }

//        xstream << "start\n";
//        for (size_t j = 0; j < branches[i]->nodes.size(); ++j)
//        {
//            auto& node0 = branches[i]->nodes[j];
//            xstream << node0->pNode.x() << " "
//                    << node0->pNode.y() << " "
//                    << node0->pNode.z() << " "
//                    << node0->radius << "\n";
//        }
//        xstream << "done\n";

//    }

//    xstream.close();


    xstream.open("/abdellah2/scratch/thinning/output/projections/radii.txt", std::ios::out);

    for (size_t i = 0; i < branches.size(); ++i)
    {
        if (branches[i]->nodes.size() == 0)
        {
            std::cout << "Empty branch \n";
            continue;
        }

        if (!branches[i]->valid)
        {
            std::cout << "Invalid branch \n";
            continue;
        }

        //xstream << "start\n";
        for (size_t j = 0; j < branches[i]->nodes.size(); ++j)
        {
            auto& node0 = branches[i]->nodes[j];
            if (node0->radius >= 2.0)
            {

                centers.push_back(Vector3f(node0->pNode.x(), node0->pNode.y(), node0->pNode.z()));
                radii.push_back(node0->radius);

                xstream << node0->pNode.x() << " "
                        << node0->pNode.y() << " "
                        << node0->pNode.z() << " "
                        << node0->radius << "\n";
            }
        }
        //xstream << "done\n";

    }

    xstream.close();



    std::vector< GraphNode* > somaRegion;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i]->radius > 1.5)
        {
            somaRegion.push_back(nodes[i]);
        }
    }



    return shellPoints;

    // Sample* somaRegion  = new Sample(somaCenter, somaRadius, 0);

    _rasterize(somaSphere, _grid);
    // return;































    // Filter the nodes
    std::vector< GraphNode* > filteredNodes;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i]->insideSoma)
            continue;

        filteredNodes.push_back(nodes[i]);
    }

    // Add the soma node to the list
    filteredNodes.push_back(somaNode);
    std::cout << "New Nodes (Samples): " << filteredNodes.size() << "\n";

    // nodes.clear();
    // nodes.shrink_to_fit();


    std::fstream sstream;
    sstream.open("/abdellah2/scratch/thinning/output/projections/file2.txt", std::ios::out);
    for (size_t i = 0; i < nodes.size(); ++i)
    {

        auto& node = nodes[i];
        sstream << node->pNode.x() << " "
               << node->pNode.y() << " "
               << node->pNode.z() << " "
               << node->radius << "\n";

    }

    sstream.close();



//    // Construct the heirarichy to the terminal
//    for (size_t i = 0; i < filteredNodes.size(); ++i)
//    {
//        //
//        auto& node = filteredNodes[i];

////        // It must be a branching node
////        if (node->connectedToSoma)
////        {
////            // Construct the branch
////            for (size_t j = 0; j < node->edgeNodes.size(); ++j)
////            {
////                auto& edgeNode = node->edgeNodes[j];

////                if (edgeNode->isSoma)
////                    continue;

////                std::cout << node->index << " " << edgeNode->index << " \n";
////                buildBranchesRecursively(node, edgeNode, branches);
////            }
////        }


//        if (node->branching)
//        {
//            // Construct the branch
//            for (size_t j = 0; j < node->edgeNodes.size(); ++j)
//            {
//                auto& edgeNode = node->edgeNodes[j];
//                {
//                    branches.push_back(buildBranch(node, edgeNode));
//                }
//            }
//        }
//    }

    std::cout << "Branches: " << branches.size() << "\n";

    std::fstream stream;
    stream.open("/abdellah2/scratch/thinning/output/projections/file.txt", std::ios::out);
    for (size_t i = 0; i < branches.size(); i++)
    {
        for (size_t j =0; j < branches[i]->nodes.size(); ++j)
        {
            auto& node = branches[i]->nodes[j];
            stream << node->pNode.x() << " "
                   << node->pNode.y() << " "
                   << node->pNode.z() << " "
                   << node->radius << "\n";
        }
    }
    stream.close();
















//    std::vector< std::vector< size_t > > branches;
//    // From every bifurcation point, start collecting the edges to construct the graph
//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        // Keep a reference to the node ID
//        const GraphNode* node = nodes[i];

//        // Only for the branching points
//        if (node->isBranching)
//        {
//            // Get all the edges connected to the node
//            for (size_t j = 0; j < node->edges.size(); ++j)
//            {
//                const GraphEdge* edge = node->edges[j];

//                std::vector< size_t > branch;
//                constructBranch(node, edge, nodes, branch);
//                // branches.push_back(branch);
//            }
//        }
//    }

//    std::cout << branches.size() << "\n : Branches \n";
}



}
