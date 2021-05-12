/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_COMMON_DEFINES_H
#define ULTRALISER_COMMON_DEFINES_H

#include <string>

#define AXES_COUNT                              6

// Refercing XYZ from arrays
#define DIMENSIONS                              3 // XYZ
#define X_DIMENSION                             0 // [0]
#define Y_DIMENSION                             1 // [1]
#define Z_DIMENSION                             2 // [2]

// Flood-filling values
#define EMPTY_VOXEL_VALUE       0
#define FILLED_VOXEL_VALUE      255

// PBRT lens radii
#define LENS_RADIUS_STEP        0.5
#define LESN_RADIUS_MAX         50

// Stream processing
#define CLEAR_STREAM stream.str(EMPTY)

// Strings
#define CHAR(C)                                 char(C)
#define STRING(STR)                             std::string(STR)
#define EMPTY                                   STRING("")
#define NO_DEFAULT_VALUE                        STRING("NO_DEFAULT_VALUE")

// Files extensions
#define HEADER_EXTENSION                        STRING(".hdr")
#define BOUNDS_EXTENSION                        STRING(".bounds")
#define RAW_EXTENSION                           STRING(".img")
#define BINARY_EXTENSION                        STRING(".bin")
#define ASCII_EXTENSION                         STRING(".ascii")
#define NRRD_EXTENSION                          STRING(".nrrd")
#define INFO_EXTENSION                          STRING(".info")
#define OMESH_INFO_EXTENSION                    STRING(".oinfo")
#define VOLUME_INFO_EXTENSION                   STRING(".vol-info")
#define MESH_INFO_EXTENSION                     STRING(".mesh-info")
#define MORPHOLOGY_INFO_EXTENSION               STRING(".morph-info")
#define WARNING_EXTENSION                       STRING(".warning")
#define VOXEL_EXTENSION                         STRING(".vox")
#define PBRT_EXTENSION                          STRING(".pbrt")
#define PNG_EXTENSION                           STRING(".png")
#define PPM_EXTENSION                           STRING(".ppm")
#define EXR_EXTENSION                           STRING(".exr")
#define OBJ_EXTENSION                           STRING(".obj")
#define OFF_EXTENSION                           STRING(".off")
#define PLY_EXTENSION                           STRING(".ply")
#define STL_EXTENSION                           STRING(".stl")

// Artifacts strings
#define INPUT_STRING                            STRING("input")
#define DMC_STRING                              STRING("dmc")
#define LAPLACIAN_STRING                        STRING("laplacian")
#define OPTIMIZED_STRING                        STRING("optimized")
#define MANIFOLD_STRING                         STRING("advanced")
#define WATERTIGHT_STRING                       STRING("watertight")
#define PROJECTION_STRING                       STEING("projection")
#define VOLUME_STRING                           STRING("volume")

// Artifacts Suffixes
#define DMC_SUFFIX                              STRING("-dmc")
#define LAPLACIAN_SUFFIX                        STRING("-laplacian")
#define OPTIMIZED_SUFFIX                        STRING("-optimized")
#define MANIFOLD_SUFFIX                         STRING("-advanced")
#define WATERTIGHT_SUFFIX                       STRING("-watertight")
#define VOLUME_MESH_SUFFIX                      STRING("-volume-mesh")
#define PROJECTION_SUFFIX                       STRING("-projection")
#define XY_SUFFIX                               STRING("-xy")
#define XZ_SUFFIX                               STRING("-xz")
#define YZ_SUFFIX                               STRING("-yz")
#define XYZ_SUFFIX                              STRING("-xyz")

// Directories where the artifacts will be generated
#define MESHES_DIRECTORY                        STRING("meshes")
#define VOLUMES_DIRECTORY                       STRING("volumes")
#define PROJECTIONS_DIRECTORY                   STRING("projections")
#define STACKS_SIRECTORY                        STRING("stacks")
#define STATISTICS_DIRECTORY                    STRING("statistics")
#define DISTRIBUTIONS_DIRECTORY                 STRING("distributions")

// Distributions stats
#define DISTRIBUTION_EXTENSION                  std::string(".dist")
#define ASPECT_RATIO_SUFFIX                     std::string("-aspect-ratio")
#define RADIUS_RATIO_SUFFIX                     std::string("-radius-ratio")
#define EDGE_RATIO_SUFFIX                       std::string("-edge-ratio")
#define RADIUS_TO_EDGE_RATIO_SUFFIX             std::string("-radius-to-edge-ratio")
#define MIN_ANGLE_SUFFIX                        std::string("-min-angle")
#define MAX_ANGLE_SUFFIX                        std::string("-max-angle")
#define TRIANGLE_SHAPE_SUFFIX                   std::string("-triangle-shape")
#define TRIANGLE_SHAPE_SIZE_SUFFIX              std::string("-triangle-size-shape")
#define SCALED_JACOBIAN_SUFFIX                  std::string("-scaled-jacobian")
#define CONDITION_NUMBER_SUFFIX                 std::string("-condition-number")
#define DISTORTION_SUFFIX                       std::string("-distortion")
#define RELATIVE_SIZE_SUFFIX                    std::string("-relative-size")
#define SAMPLES_RADII                           std::string("-samples-radii")
#define SECTION_AVERAGE_RADIUS                  std::string("-section-average-radius")
#define NUMBER_SAMPLES_PER_SECTION              std::string("-number-samples-per-section")
#define SEGMENTS_LENGTH                         std::string("-segments-length")
#define SECTIONS_LENGTH                         std::string("-sections-length")
#define SEGMENTS_SURFACE_AREA                   std::string("-segments-surface-area")
#define SECTIONS_SURFACE_AREA                   std::string("-sections-surface-area")
#define SEGMENTS_VOLUME                         std::string("-segments-volume")
#define SECTIONS_VOLUME                         std::string("-sections-volume")


// OBJ flags
#define OBJ_VERTEX_FLAG                         STRING("v")
#define OBJ_VERTEX_NORMAL_FLAG                  STRING("vn")
#define OBJ_FACE_FLAG                           STRING("f")
#define OBJ_TEXTURE_FLAG                        STRING("vt")
#define OBJ_HASH                                CHAR('#')
#define OBJ_END_FLAG                            STRING("#end")

// PLY flags
#define PLY_VERTEX_FLAG                         STRING("element vertex")
#define PLY_FACE_FLAG                           STRING("element face")
#define PLY_TEXTURE_FLAG                        STRING("property float s")
#define PLY_END_FLAG                            STRING("end_header")
#define PLY_HEADER                              STRING("ply")
#define PLY_HSIZE                               3
#define PLY_FORMAT_ASCII                        0
#define PLY_FORMAT_BIN_L                        1
#define PLY_FORMAT_BIN_B                        2

// STL glags
#define STL_SOLID_KEYWORD                       STRING("solid")
#define STL_FACET_KEYWORD                       STRING("facet")
#define STL_VERTEX_KEYWORD                      STRING("vertex")

// OFF flags
#define OFF_HEADER                              STRING("OFF")
#define OFF_HSIZE                               3


// String flags
#define BACK_SLASH                              STRING("/")
#define DOUBLE_BACK_SLASH                       STRING("//")
#define SPACE                                   STRING(" ")
#define NEW_LINE                                STRING("\n")
#define C_SPACE                                 CHAR(' ')
#define C_BACK_SLASH                            CHAR('/')

// Output directories
#define DOF_DIRECTORY                           STRING("dof")
#define XY_DIRECTORY                            STRING("xy")
#define ZY_DIRECTORY                            STRING("zy")
#define BRIGHTFIELD_DIRECTORY                   STRING("pbrt-brightfield")
#define FLUORESCENCE_DIRECTORY                  STRING("pbrt-fluorescence")
#define LSFM_DIRECTORY                          STRING("pbrt-lsfm")
#define FOCUS_ALL_OBJECTS                       STRING("all-objects")
#define FOCUS_ALL_SLICES                        STRING("all-slices")
#define FIXED_FOCUS_LR                          STRING("fixed-focus-lr")

// Configuration keys (for key-value mapping)
#define GID_KEY                                 STRING("GID:")
#define MTYPE_KEY                               STRING("MORPHOLOGY_TYPE:")
#define MLABEL_KEY                              STRING("MORPHOLOGY_LABEL:")
#define TAG_KEY                                 STRING("TAG:")
#define POSITION_KEY                            STRING("POSITION:")
#define ORIENTATION_KEY                         STRING("ORIENTATION:")
#define TRANSFORM_KEY                           STRING("TRANSFORM:")
#define COLUMN_KEY                              STRING("COLUMN:")
#define LAYER_KEY                               STRING("LAYER:")

#define AXON_INDEX                              1
#define APICAL_DENDRITE_INDEX                   2
#define BASAL_DENDRITE_INDEX                    3
#define SOMA_INDEX                              4

// Logging
#define TIME_STAMP_CHAR_LENGTH                  128
#define PROGRESS_BAR_LENGTH                     50
#define TITLE_LENGTH                            80

// Application done
#define ULTRALISER_DONE                         printf("\n"); return EXIT_SUCCESS;

#endif // ULTRALISER_COMMON_DEFINES_H
