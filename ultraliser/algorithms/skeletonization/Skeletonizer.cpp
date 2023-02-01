#include "Skeletonizer.h"
#include <math/Vector.h>

namespace Ultraliser
{
Skeletonizer::Skeletonizer(const Mesh *mesh,
                           Volume* volume)
    : _mesh(mesh)
    , _volume(volume)
{
    // Mesh bounding box
    _pMinMesh = volume->getPMin();
    _pMaxMesh = volume->getPMax();
    _boundsMesh = _pMaxMesh - _pMinMesh;
    _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

    // Volume bounding box
    _pMinVolume = Vector3f(0.f);
    _pMaxVolume = Vector3f(volume->getWidth() * 1.f,
                           volume->getHeight() * 1.f,
                           volume->getDepth() * 1.f);
    _boundsVolume = _pMaxVolume;
    _centerVolume = 0.5 * _boundsVolume;

    // Mesh to volume scale factor
    _scaleFactor = _boundsMesh / _boundsVolume;
}

void Skeletonizer::applyVolumeThinning()
{
    TIMER_SET;
    LOG_STATUS("Thinning Volume");

    // A block of the volume that is scanned every iteration
    int8_t volumeBlock[26];

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel = std::make_unique<Thinning6Iterations>();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    LOOP_STARTS("Thinning Loop");
    LOOP_PROGRESS(0, 100);
    while(1)
    {
        size_t numberDeletedVoxels = 0;

        // Search for the border voxels, for this iteration
        std::vector< std::vector< Vec3ui_64 > > borderVoxels = _volume->searchForBorderVoxels();

        for (size_t direction = 0; direction < 6; direction++)
        {
            // Search for the delerable voxels
            std::vector< Vec3ui_64 > voxelsToBeDeleted =
                    _volume->searchForDeletableVoxels(
                        borderVoxels, thinningKernel, direction, volumeBlock);

            // Delete the voxels
            for (size_t i = 0; i < voxelsToBeDeleted.size(); ++i)
            {
                numberDeletedVoxels++;
                _volume->clear(voxelsToBeDeleted[i].x(),
                               voxelsToBeDeleted[i].y(),
                               voxelsToBeDeleted[i].z());
            }

            // Clear the container of the deleted voxels
            voxelsToBeDeleted.clear();
        }

        // Clear the border voxels (list of lists)
        for (size_t i = 0; i < borderVoxels.size(); ++i)
        {
            borderVoxels[i].clear();
            borderVoxels[i].shrink_to_fit();
        }
        borderVoxels.clear();
        borderVoxels.shrink_to_fit();

         // Updating the progess bar
        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
        LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                      initialNumberVoxelsToBeDeleted);

        if (numberDeletedVoxels == 0)
            break;

        loopCounter++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

std::vector< Vector3f > Skeletonizer::getShellPoints()
{
    if (_shellPoints.size() == 0)
    {
        _computeShellPoints();
    }

    return _shellPoints;
}

void Skeletonizer::_computeShellPoints()
{
    // Search for the border voxels (the shell voxels) of the volume
    std::vector< std::vector< Vec3ui_64 > > perSliceSurfaceShell = _volume->searchForBorderVoxels();

    // Concatinate the points in a single list
    for (size_t i = 0; i < perSliceSurfaceShell.size(); ++i)
    {
        for (size_t j = 0; j < perSliceSurfaceShell[i].size(); ++j)
        {
            const auto voxel = perSliceSurfaceShell[i][j];
            _shellPoints.push_back(Vector3f(voxel.x(), voxel.y(), voxel.z()));
        }
        perSliceSurfaceShell[i].clear();
        perSliceSurfaceShell[i].shrink_to_fit();
    }
    perSliceSurfaceShell.clear();
    perSliceSurfaceShell.shrink_to_fit();

    // Adjust the locations of the shell points taking into consideration the mesh coordinates
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _shellPoints.size(); ++i)
    {
        // Center the shell points (of the volume) at the origin
        _shellPoints[i] -= _centerVolume;

        // Scale to match the dimensions of the mesh
        _shellPoints[i].x() *= _scaleFactor.x();
        _shellPoints[i].y() *= _scaleFactor.y();
        _shellPoints[i].z() *= _scaleFactor.z();

        // Translate to the center of the mesh
        _shellPoints[i] += _centerMesh;
    }
}



}
