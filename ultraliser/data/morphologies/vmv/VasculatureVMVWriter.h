#ifndef VASCULATUREVMVWRITER_H
#define VASCULATUREVMVWRITER_H

#include <data/morphologies/VasculatureMorphology.h>

namespace Ultraliser
{

/**
 * @brief writeMorphologyToVMVFile
 * @param morphology
 * @param prefix
 */
void writeMorphologyToVMVFile(const VasculatureMorphology *morphology, const std::string &prefix);

}


#endif // VASCULATUREVMVWRITER_H
