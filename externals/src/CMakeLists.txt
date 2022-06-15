# Ultraliser library
set(NRRDREADER_LIBRARY "nrrdreader")

# Collect the header and source files
file(GLOB_RECURSE NRRD_READER_HEADERS "*.h")
file(GLOB_RECURSE NRRD_READER_SOURCES "*.cpp")

# Add the current directory to include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Create the library
add_library(${NRRDREADER_LIBRARY} ${NRRD_READER_HEADERS} ${NRRD_READER_SOURCES})

# Copy the artifacts for the apps
set_property(TARGET ${NRRDREADER_LIBRARY} PROPERTY IMPORTED_LOCATION
    ${CMAKE_BINARY_DIRECTORY}/lib${NRRDREADER_LIBRARY}.so)

# Install the files
install(TARGETS ${NRRDREADER_LIBRARY} DESTINATION lib)
install(FILES ${NRRDREADER_LIBRARY} DESTINATION include)

# Link the libraries
target_link_libraries(${NRRDREADER_LIBRARY} PRIVATE ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})
target_link_libraries(${NRRDREADER_LIBRARY} PUBLIC fmt)

target_include_directories(nrrdreader PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
    "$<INSTALL_INTERFACE:include>")

target_compile_options(nrrdreader PRIVATE -Wall -Wextra -pedantic -Werror)
