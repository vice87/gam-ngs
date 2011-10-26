# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/api/
# ==========================

# list include paths
include_directories( ${BamTools_SOURCE_DIR}/src )

# add compiler definitions 
add_definitions( -DBAMTOOLS_API_LIBRARY ) # (for proper exporting of library symbols)
add_definitions( -fPIC ) # (attempt to force PIC compiling on CentOS, not being set on shared libs by CMake)

# list of all BamTools API source (.cpp) files
set( BamToolsAPISources
        BamAlignment.cpp
        BamMultiReader.cpp
        BamReader.cpp
        BamWriter.cpp
        SamHeader.cpp
        SamProgram.cpp
        SamProgramChain.cpp
        SamReadGroup.cpp
        SamReadGroupDictionary.cpp
        SamSequence.cpp
        SamSequenceDictionary.cpp
        internal/BamDeviceFactory_p.cpp
        internal/BamException_p.cpp
        internal/BamFile_p.cpp
        internal/BamFtp_p.cpp
        internal/BamHeader_p.cpp
        internal/BamHttp_p.cpp
        internal/BamIndexFactory_p.cpp
        internal/BamMultiReader_p.cpp
        internal/BamPipe_p.cpp
        internal/BamRandomAccessController_p.cpp
        internal/BamReader_p.cpp
        internal/BamStandardIndex_p.cpp
        internal/BamToolsIndex_p.cpp
        internal/BamWriter_p.cpp
        internal/BgzfStream_p.cpp
        internal/ILocalIODevice_p.cpp
        internal/IRemoteIODevice_p.cpp
        internal/SamFormatParser_p.cpp
        internal/SamFormatPrinter_p.cpp
        internal/SamHeaderValidator_p.cpp
)

set (CMAKE_BUILD_TYPE Release)

# create BamTools static library
add_library( BamTools STATIC ${BamToolsAPISources} )
target_link_libraries( BamTools ${ZLIB_LIBRARIES})

set_target_properties( BamTools PROPERTIES OUTPUT_NAME "bamtools" )
set_target_properties( BamTools PROPERTIES PREFIX "lib" )
