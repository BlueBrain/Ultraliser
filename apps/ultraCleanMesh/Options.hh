#ifndef OPTIONS_HH
#define OPTIONS_HH

#include <common/Common.h>

/**
 * @brief The Options struct
 */
struct Options
{
    /**
     * @brief inputMesh
     * Input mesh file
     */
    std::string inputMesh;

    /**
     * @brief outputDirectory
     * The directory where the volume will be created.
     * Output directory
     */
    std::string outputDirectory;

    /**
     * @brief boundsFile
     * Use a bounds file to only voxelize part of the mesh or even a greater
     * space. If the file is not given, the bounding box of the input mesh will
     * be used in addition to a little delta to avoid intersection.
     */
    std::string boundsFile;

    /**
     * @brief volumeResolution
     * The base resolution of the volume that corresponds to the largest
     * dimension.
     */
    uint64_t volumeResolution;

    /**
     * @brief autoResolution
     * Sets the resolution of the volume based on mesh dimensions.
     */
    bool autoResolution;

    /**
     * @brief voxelsPerMicron
     * Number of voxels per micron in case of auto resolution.
     */
    uint64_t voxelsPerMicron;

    /**
     * @brief edgeGap
     */
    float edgeGap;

    /**
     * @brief solid
     * Fill the interior of the volume using solid voxelization.
     */
    bool solid;

    /**
     * @brief volumeType
     * Use a specific volume for the voxelization process.
     */
    std::string volumeType;

    /**
     * @brief projectXY
     * If this flag is set, the XY projection of the volume will be saved to a
     * PNG image. This flag is set to validate the output volume.
     */
    bool projectXY;

    /**
     * @brief projectZY
     * If this flag is set, the ZY projection of the volume will be saved to a
     * PNG image. This flag is set to validate the output volume.
     */
    bool projectZY;

    /**
     * @brief stackXY
     * Create an image stack along the XY plane.
     */
    bool stackXY;

    /**
     * @brief stackZY
     * Create an image stack along the ZY plane.
     */
    bool stackZY;

    /**
     * @brief createBinaryVolume
     * If this flag is set, a binary volume will be created. This volume has
     * 1 bit per voxel.
     */
    bool writeBitVolume;

    /**
     * @brief createByteVolume
     * If this flag is set, a default raw volume will be created. This volume
     * has 1 byte per voxel.
     */
    bool writeByteVolume;

    /**
     * @brief writeNRRDVolume
     * If this flag is set, the volume will be written to an NRRD file that is
     * compatible with VTK.
     */
    bool writeNRRDVolume;

    /**
     * @brief voxelizeMesh
     * If this flag is set, the input mesh will be voxelized. If any of the
     * volume export flags are set, the volume resulting from this process will
     * be written to disk.
     */
    bool voxelizeMesh;

    /**
     * @brief optimizationIterations
     * Number of iterations of optimizing the mesh surface.
     * By default it is set to 1, but more accurate numbers depend on the
     * given mesh with trial and error.
     */
    float optimizationIterations;

    /**
     * @brief smoothingIterations
     * Number of iterations required to smooth the optimized mesh, by default 10.
     */
    int64_t smoothingIterations;

    /**
     * @brief flatFactor
     * A factor that is used for the coarseFlat function.
     * Default value is 0.1.
     */
    float flatFactor;

    /**
     * @brief denseFactor
     * A factor that is used for the coarseDense function.
     * Default value is 5.0.
     */
    float denseFactor;

    /**
     * @brief exportOBJ
     * Export any reconstructed mesh to .OBJ file.
     */
    bool exportOBJ;

    /**
     * @brief exportPLY
     * Export any reconstructed mesh to .PLY file.
     */
    bool exportPLY;

    /**
     * @brief exportOFF
     * Export any reconstructed mesh to .OFF file.
     */
    bool exportOFF;

    /**
     * @brief exportSTL
     * Export any reconstructed mesh to .STL file.
     */
    bool exportSTL;

    /**
     * @brief prefix
     * Just a prefix that will be used to label the output files. If this
     * is not given by the user, the name of the mesh file will be used.
     */
    std::string prefix;

    /**
     * @brief writeStatictics
     * Write the statictics.
     */
    bool writeStatistics;

    /**
     * @brief outputPrefix
     * Simply, the [OUTPUT_DIRECTORY]/[PREFIX]. This variable is just added to
     * make the code simpler.
     */
    std::string outputPrefix;
};

#endif // OPTIONS_HH