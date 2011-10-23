cmake_minimum_required (VERSION 2.6)
project(GamNgs)

add_definitions( -Wno-deprecated )

find_package(Boost COMPONENTS graph REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories("${PROJECT_SOURCE_DIR}/lib/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/BamFile/include")
include_directories( ${Boost_INCLUDE_DIRS} )

# sorgenti da compilare
file(GLOB GAMNGSLIB_SRC_FILES 
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/ablast.cc
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/aligner.cc
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/alignment.cc
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/smith_waterman.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/contig.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/io_contig.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/nucleotide.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/Read.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/Frame.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/Block.cc
    ${PROJECT_SOURCE_DIR}/lib/src/graphs/AssemblyGraph.cc
    ${PROJECT_SOURCE_DIR}/lib/src/graphs/PairingEvidencesGraph.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/BestPctgCtgAlignment.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/ContigInPctgInfo.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/PairedContig.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/PctgBuilder.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pool/HashContigMemPool.cc
    ${PROJECT_SOURCE_DIR}/lib/src/strand_fixer/RelativeStrand.cc
    ${PROJECT_SOURCE_DIR}/lib/src/strand_fixer/RelativeStrandEvidences.cc
    ${PROJECT_SOURCE_DIR}/lib/src/strand_fixer/StrandProbability.cc
    ${PROJECT_SOURCE_DIR}/lib/src/OrderingFunctions.cc
    ${PROJECT_SOURCE_DIR}/lib/src/PartitionFunctions.cc
)

# libreria contenente le classi usate
#add_subdirectory(lib)

# add the executable
add_executable(gam-ngs src/gam-ngs.cc ${GAMNGSLIB_SRC_FILES})

target_link_libraries(gam-ngs z)
target_link_libraries(gam-ngs ${PROJECT_SOURCE_DIR}/lib/BamFile/libbamtools.a)
target_link_libraries(gam-ngs ${Boost_LIBRARIES})

# executable for computing fasta statistics
add_executable(n50 src/n50.cc)

# ----------------------------------------------------------
# - FOLLOWING EXECUTABLES ARE FOR TESTING ONLY
# ----------------------------------------------------------

# executable for memory usage debugging
#add_executable(memory-usage src/test/memory-usage.cc 
#    lib/src/pool/HashContigMemPool.cc 
#    lib/src/pool/ContigMemPool.cc 
#    lib/src/assembly/contig.cc 
#    lib/src/assembly/io_contig.cc 
#    lib/src/assembly/nucleotide.cc
#    lib/src/alignment/alignment.cc
#)
#target_link_libraries(memory-usage z)
#target_link_libraries(memory-usage ${PROJECT_SOURCE_DIR}/lib/BamFile/libbamtools.a)

# executable for debugging alignment in weaving
#add_executable(alignment-debug src/test/alignment.cc ${GAMNGSLIB_SRC_FILES})
#target_link_libraries(alignment-debug z)
#target_link_libraries(alignment-debug ${PROJECT_SOURCE_DIR}/lib/BamFile/libbamtools.a)

#add_executable(filter-blocks src/filter-blocks.cc ${GAMNGSLIB_SRC_FILES})
#target_link_libraries(filter-blocks "${PROJECT_SOURCE_DIR}/BamFile/lib/libbamtools.so")
#target_link_libraries(filter-blocks  ${Boost_LIBRARIES})