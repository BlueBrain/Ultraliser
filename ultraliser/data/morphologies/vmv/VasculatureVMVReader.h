#ifndef VASCULATURE_VMV_READER_H
#define VASCULATURE_VMV_READER_H

#include <data/morphologies/h5/VasculatureH5Sample.hh>
#include <data/morphologies/h5/VasculatureH5Section.hh>
#include <data/morphologies/h5/VasculatureH5Connectivity.hh>
#include <data/morphologies/VasculatureMorphology.h>

namespace Ultraliser
{

class VasculatureVMVReader
{
public:

    /**
     * @brief VasculatureVMVReader
     * Constructor
     * @param vmvMorphologyFilePath
     * The path to the VMV morphology file.
     */
    VasculatureVMVReader(const std::string &vmvMorphologyFilePath);

public:

    /**
     * @brief getMorphology
     * Return a pointer to the vasculature morphology.
     * @return Return a pointer to the morphology.
     */
    VasculatureMorphology* getMorphology();

private:

    /**
     * @brief _readAttributes
     * Reads the atttributes from the morphology file to be able to process the file easily.
     */
    void _readAttributes();

    /**
     * @brief _readSamples
     * Reads the samples from the morphology file.
     */
    void _readSamples();

    /**
     * @brief _readStrands
     * Reads the connectivity information form the morphology file.
     */
    void _readStrands();

private:

    /**
     * @brief _numberVerts
     * Number of samples in the morphology.
     */
    uint64_t _numberVerts;

    /**
     * @brief _numberStrands
     * Number of strands or edges in the morphology.
     */
    uint64_t _numberStrands;

    /**
     * @brief _numberAttributesPerVertex
     * Number of attributes per vertex.
     */
    uint64_t _numberAttributesPerVertex;

    /**
     * @brief _vmvMorphologyFile
     * The path to the VMV morphology file.
     */
    std::string _vmvMorphologyFile;

    /**
     * @brief _samples
     * Vasculature samples.
     */
    Samples _samples;

    /**
     * @brief _sections
     * Vasculature sections or strands.
     */
    Sections _sections;
};

}

#endif // VASCULATURE_VMV_READER_H
