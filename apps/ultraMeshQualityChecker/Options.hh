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
     * @brief prefix
     * Just a prefix that will be used to label the output files. If this
     * is not given by the user, the name of the mesh file will be used.
     */
    std::string prefix;

    /**
     * @brief outputPrefix
     * Simply, the [OUTPUT_DIRECTORY]/[PREFIX]. This variable is just added to
     * make the code simpler.
     */
    std::string outputPrefix;
};

#endif // OPTIONS_HH
