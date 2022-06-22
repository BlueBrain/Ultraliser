#include "Arrays.h"
#include <cstddef>
#include <cstring>


namespace Ultraliser
{
namespace Array
{

uint8_t* convertStringTo8UIArray(const std::string &stringArray)
{
    // Array size
    const auto arraySize = stringArray.size();

    // Create the data array
    uint8_t* data = new uint8_t[arraySize];

    // Fill the data array
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < arraySize; ++i)
    {
        data[i] = stringArray[i];
    }

    // Return the data array
    return data;
}

uint16_t* convertStringTo16UIArray(const std::string &stringArray)
{
    // Array size
    const auto arraySize = static_cast< size_t >(stringArray.size() / 2.0);

    // Create the data array
    uint16_t* data = new uint16_t[arraySize];

    // Fill the data array
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < arraySize; ++i)
    {
        uint8_t v0 = stringArray[2 * i];
        uint8_t v1 = stringArray[2 * i + 1];

        uint16_t value = 0;
        value |= static_cast< uint16_t >(v0) << 0;
        value |= static_cast< uint16_t >(v1) << 8;
        data[i] = value;
    }

    // Return the data array
    return data;
}

uint32_t* convertStringTo32UIArray(const std::string &stringArray)
{
    // Array size
    const auto arraySize = static_cast< uint64_t >(stringArray.size() / 4.0);

    // Create the data array
    uint32_t* data = new uint32_t[arraySize];

    // Fill the data array
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < arraySize; ++i)
    {
        uint8_t v0 = stringArray[4 * i];
        uint8_t v1 = stringArray[4 * i + 1];
        uint8_t v2 = stringArray[4 * i + 2];
        uint8_t v3 = stringArray[4 * i + 3];

        uint32_t value = 0;
        value |= static_cast< uint32_t >(v0) << 0;
        value |= static_cast< uint32_t >(v1) << 8;
        value |= static_cast< uint32_t >(v2) << 16;
        value |= static_cast< uint32_t >(v3) << 24;
        data[i] = value;
    }

    // Return the data array
    return data;
}

uint64_t* convertStringTo64UIArray(const std::string &stringArray)
{
    // Array size
    const auto arraySize = static_cast< uint64_t >(stringArray.size() / 8.0);

    // Create the data array
    uint64_t* data = new uint64_t[arraySize];

    // Fill the data array
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < arraySize; ++i)
    {
        uint8_t v0 = stringArray[8 * i];
        uint8_t v1 = stringArray[8 * i + 1];
        uint8_t v2 = stringArray[8 * i + 2];
        uint8_t v3 = stringArray[8 * i + 3];
        uint8_t v4 = stringArray[8 * i + 4];
        uint8_t v5 = stringArray[8 * i + 5];
        uint8_t v6 = stringArray[8 * i + 6];
        uint8_t v7 = stringArray[8 * i + 7];

        uint64_t value = 0;
        value |= static_cast< uint64_t >(v0) << 0;
        value |= static_cast< uint64_t >(v1) << 8;
        value |= static_cast< uint64_t >(v2) << 16;
        value |= static_cast< uint64_t >(v3) << 24;
        value |= static_cast< uint64_t >(v4) << 32;
        value |= static_cast< uint64_t >(v5) << 40;
        value |= static_cast< uint64_t >(v6) << 48;
        value |= static_cast< uint64_t >(v7) << 56;
        data[i] = value;
    }

    // Return the data array
    return data;
}

BitArray* convertStringToBitArray(const std::string &stringArray, const size_t &numberBits)
{
    // Initially construct the data array
    uint8_t* dataArray = convertStringTo8UIArray(stringArray);

    // Construct the bit array from the binary data array
    BitArray* bitArray = new BitArray(dataArray, numberBits);

    // Release the data array
    delete [] dataArray;

    // Return the bit array
    return bitArray;
}

}
}
