#pragma once

#include "Volume.h"
#include <utilities/Timer.h>
#include <data/volumes/grids/Grids.h>

template
void writeUnsignedGridToNRRDFile< uint8_t >(const std::string &, const VolumeGridU8*);

template
void writeUnsignedGridToNRRDFile< uint16_t >(const std::string &, const VolumeGridU16*);

template
void writeUnsignedGridToNRRDFile< uint32_t >(const std::string &, const VolumeGridU32*);

template
void writeUnsignedGridToNRRDFile< uint64_t >(const std::string &, const VolumeGridU64*);

template
void writeFloatGridToNRRDFile< float >(const std::string&, const VolumeGridF32*);

template
void writeFloatGridToNRRDFile< double >(const std::string&, const VolumeGridF64*);

template
void writeUnsignedGridToVOLFile< uint8_t >(const std::string &, const VolumeGridU8*);

template
void writeUnsignedGridToVOLFile< uint16_t >(const std::string &, const VolumeGridU16*);

template
void writeUnsignedGridToVOLFile< uint32_t >(const std::string &, const VolumeGridU32*);

template
void writeUnsignedGridToVOLFile< uint64_t >(const std::string &, const VolumeGridU64*);

template
void writeFloatGridToVOLFile< float >(const std::string&, const VolumeGridF32*);

template
void writeFloatGridToVOLFile< double >(const std::string&, const VolumeGridF64*);
