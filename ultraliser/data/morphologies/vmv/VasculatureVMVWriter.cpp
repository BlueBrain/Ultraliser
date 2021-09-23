#include "VasculatureVMVWriter.h"
#include <common/Common.h>

namespace Ultraliser
{

void writeMorphologyToVMVFile(const VasculatureMorphology* morphology, const std::string& prefix)
{
    // File
    std::ofstream stream;
    const std::string path = prefix + VMV_EXTENSION;
    stream.open(path.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot open the file [ %s ]", path.c_str());
    }

    // A reference to the samples
    auto samples = morphology->getSamples();
    auto sections = morphology->getSections();

    // Write the header
    stream << "$PARAM_BEGIN" << std::endl;
    stream << "$NUM_VERTS " << samples.size() << std::endl;
    stream << "$NUM_STRANDS " << sections.size() << std::endl;
    stream << "$PARAM_END" << std::endl << std::endl;

    // Write the samples
    stream << "$VERT_LIST_BEGIN" << std::endl;
    for (const auto& sample: samples)
    {
        const auto& position = sample->getPosition();
        const auto& radius = sample->getRadius();
        stream << sample->getIndex() << " "
               << position.x() << " " << position.y() << " " << position.z()
               << " " << radius << std::endl;
    }
    stream << "$VERT_LIST_END" << std::endl << std::endl;

    // Write the strands
    stream << "$STRANDS_LIST_BEGIN" << std::endl;
    for (const auto& section: sections)
    {
        stream << section->getIndex() << " ";
        for (const auto& samples: section->getSamples())
        {
            stream << samples->getIndex() << " ";
        }
        stream << std::endl;
    }
    stream << "$STRANDS_LIST_END" << std::endl << std::endl;

    // Write the connectivity
    stream << "$CONNECTIVITY_LIST_BEGIN" << std::endl;
    for (const auto& section: sections)
    {
        stream << section->getIndex() << " ";
        if (section->getParentIndices().size() == 0)
        {
            stream << "None ";
        }
        else
        {
            for (const auto& parent: section->getParentIndices())
            {
                stream << parent + 1 << " ";
            }
        }

        // Add a separating /
        stream << "/";

        if (section->getChildrenIndices().size() == 0)
        {
            stream << " None" << std::endl;
        }
        else
        {
            for (const auto& child: section->getChildrenIndices())
            {
                stream << " " << child + 1;
            }
            stream << std::endl;
        }
    }
    stream << "$CONNECTIVITY_LIST_END" << std::endl << std::endl;

    // Close the file
    stream.close();
}

}
