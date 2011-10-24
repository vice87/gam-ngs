# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/api/
# ==========================

# list include paths
# include_directories( ${BamTools_SOURCE_DIR}/src )

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
        internal/BamHeader_p.cpp
        internal/BamIndexFactory_p.cpp
        internal/BamMultiReader_p.cpp
        internal/BamRandomAccessController_p.cpp
        internal/BamReader_p.cpp
        internal/BamStandardIndex_p.cpp
        internal/BamToolsIndex_p.cpp
        internal/BamWriter_p.cpp
        internal/BgzfStream_p.cpp
        internal/SamFormatParser_p.cpp
        internal/SamFormatPrinter_p.cpp
        internal/SamHeaderValidator_p.cpp
)

# create BamTools static library
add_library( BamTools STATIC ${BamToolsAPISources} )
target_link_libraries( BamTools ${ZLIB_LIBRARIES})
set_target_properties( BamTools PROPERTIES OUTPUT_NAME "bamtools" )
set_target_properties( BamTools PROPERTIES PREFIX "lib" )