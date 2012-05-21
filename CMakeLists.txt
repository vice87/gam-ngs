cmake_minimum_required (VERSION 2.6)
project(GamNgs)

add_definitions( -Wno-deprecated )

find_package(Boost COMPONENTS graph program_options system filesystem REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories("${PROJECT_SOURCE_DIR}/lib/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/bamtools-2.0.5/src")
include_directories("${PROJECT_SOURCE_DIR}/lib/murmur3")
include_directories( ${Boost_INCLUDE_DIRS} )

# sorgenti da compilare
file(GLOB GAMNGSLIB_SRC_FILES
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/ablast_new.cc
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/my_alignment.cc
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/full_smith_waterman.cc
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/banded_smith_waterman.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/contig.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/io_contig.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/nucleotide.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/Read.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/Frame.cc
    ${PROJECT_SOURCE_DIR}/lib/src/assembly/Block.cc
    ${PROJECT_SOURCE_DIR}/lib/src/bam/MultiBamReader.cc
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
    ${PROJECT_SOURCE_DIR}/lib/src/PartitionFunctions.cc
    ${PROJECT_SOURCE_DIR}/lib/src/UtilityFunctions.cc
    ${PROJECT_SOURCE_DIR}/lib/murmur3/murmur3.cc
)

# sottocartella dove si trova la libreria BamTools
add_subdirectory(lib)

# GAM-NGS EXECUTABLE
add_executable(gam-ngs src/gam-exp.cc src/Merge.cc src/CreateBlocks.cc src/Options.cc src/OptionsStandard.cc ${GAMNGSLIB_SRC_FILES})

target_link_libraries(gam-ngs ${ZLIB_LIBRARIES})
target_link_libraries(gam-ngs ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gam-ngs BamTools) #target_link_libraries(gam-ngs ${PROJECT_SOURCE_DIR}/lib/BamFile/libbamtools.a)
target_link_libraries(gam-ngs ${Boost_LIBRARIES})

# GAM-NGS FOR FOSMID POOLS EXECUTABLE
#add_executable(gam-fp src/gam-fp.cc ${GAMNGSLIB_SRC_FILES})
#
#target_link_libraries(gam-fp ${ZLIB_LIBRARIES})
#target_link_libraries(gam-fp ${CMAKE_THREAD_LIBS_INIT})
#target_link_libraries(gam-fp BamTools)
#target_link_libraries(gam-fp ${Boost_LIBRARIES})

# executable for computing fasta statistics
add_executable(n50 src/n50.cc)

# ----------------------------------------------------------
# - FOLLOWING EXECUTABLES ARE FOR TESTING ONLY
# ----------------------------------------------------------

# executable to get number of paired-end reads correctly placed.

# executable for getting unused contigs from "old-format" pctgs file
add_executable(pe-stats src/pe-stats.cc)
target_link_libraries(pe-stats BamTools)

add_executable(hash-test src/test/hash_test.cc lib/src/bam/MultiBamReader.cc lib/src/assembly/Read.cc lib/src/UtilityFunctions.cc lib/murmur3/murmur3.cc)
target_link_libraries(hash-test BamTools)

add_executable(pe-test src/test/pe-test.cc)
target_link_libraries(pe-test BamTools)

add_executable(bam-stats lib/src/UtilityFunctions.cc lib/src/bam/MultiBamReader.cc src/test/bam-stats.cc)
target_link_libraries(bam-stats BamTools)

#add_executable(blocks-coverage lib/src/UtilityFunctions.cc lib/src/bam/MultiBamReader.cc lib/src/assembly/Block.cc lib/src/assembly/Frame.cc lib/src/assembly/Read.cc src/test/blocks-coverage.cc)
#target_link_libraries(blocks-coverage BamTools)

# executable for testing new aligners
#add_executable(new-aligners src/test/new-alignments.cc
#    lib/src/assembly/contig.cc
#    lib/src/assembly/io_contig.cc
#    lib/src/assembly/nucleotide.cc
#    lib/src/alignment/my_alignment.cc
#    lib/src/alignment/full_smith_waterman.cc
#    lib/src/alignment/banded_smith_waterman.cc
#    lib/src/alignment/ablast_new.cc
#)

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
