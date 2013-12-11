#  This file is part of GAM-NGS.
#  Copyright (c) 2011 by Riccardo Vicedomini <rvicedomini@appliedgenomics.org>,
#  Francesco Vezzi <vezzi@appliedgenomics.org>,
#  Simone Scalabrin <scalabrin@appliedgenomics.org>,
#  Lars Arverstad <lars.arvestad@scilifelab.se>,
#  Alberto Policriti <policriti@appliedgenomics.org>,
#  Alberto Casagrande <casagrande@appliedgenomics.org>
#
#  GAM-NGS is an evolution of a previous work (GAM) done by Alberto Casagrande,
#  Cristian Del Fabbro, Simone Scalabrin, and Alberto Policriti.
#  In particular, GAM-NGS has been adapted to work on NGS data sets and it has
#  been written using GAM's software as starting point. Thus, it shares part of
#  GAM's source code.
#
#  Moreover, GAM-NGS uses BamTools library to access BAM files.
#  BamTools's source code has been put in ./lib/bamtools-2.0.5/ folder,
#  in which its license can be found.
#
#  GAM-NGS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GAM-NGS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with GAM-NGS.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.6)
project(GamNgs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_definitions( -Wno-deprecated )

find_package(Boost COMPONENTS graph program_options system filesystem REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(Sparsehash REQUIRED)

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )

# define compiler flags for all code
set( CMAKE_BUILD_TYPE Release )

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories("${PROJECT_SOURCE_DIR}/lib/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/bamtools-2.0.5/src")
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${SPARSEHASH_INCLUDE_DIR} )

# sorgenti da compilare
file(GLOB GAMNGSLIB_SRC_FILES
    ${PROJECT_SOURCE_DIR}/lib/src/alignment/ablast.cc
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
	${PROJECT_SOURCE_DIR}/lib/src/graphs/CompactAssemblyGraph.cc
    ${PROJECT_SOURCE_DIR}/lib/src/graphs/PairingEvidencesGraph.cc
	${PROJECT_SOURCE_DIR}/lib/src/pctg/BestCtgAlignment.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/BestPctgCtgAlignment.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/ContigInPctgInfo.cc
	${PROJECT_SOURCE_DIR}/lib/src/pctg/CtgInPctgInfo.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/PairedContig.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pctg/PctgBuilder.cc
	${PROJECT_SOURCE_DIR}/lib/src/pctg/ThreadedBuildPctg.cc
	${PROJECT_SOURCE_DIR}/lib/src/pctg/BuildPctgFunctions.cc
    ${PROJECT_SOURCE_DIR}/lib/src/pool/HashContigMemPool.cc
    ${PROJECT_SOURCE_DIR}/lib/src/strand_fixer/RelativeStrand.cc
    ${PROJECT_SOURCE_DIR}/lib/src/strand_fixer/RelativeStrandEvidences.cc
    ${PROJECT_SOURCE_DIR}/lib/src/strand_fixer/StrandProbability.cc
    ${PROJECT_SOURCE_DIR}/lib/src/PartitionFunctions.cc
    ${PROJECT_SOURCE_DIR}/lib/src/UtilityFunctions.cc
)

# sorgenti da compilare
file(GLOB GAM_CREATE_LIB_SRC_FILES
	${PROJECT_SOURCE_DIR}/lib/src/assembly/Read.cc
	${PROJECT_SOURCE_DIR}/lib/src/assembly/Frame.cc
	${PROJECT_SOURCE_DIR}/lib/src/assembly/Block.cc
	${PROJECT_SOURCE_DIR}/lib/src/bam/MultiBamReader.cc
	${PROJECT_SOURCE_DIR}/lib/src/UtilityFunctions.cc
)

# sottocartella dove si trova la libreria BamTools
add_subdirectory(lib)

# GAM-CREATE executable
add_executable(gam-create src/gam-create.cc src/CreateBlocks.cc src/Options.cc src/OptionsCreate.cc ${GAM_CREATE_LIB_SRC_FILES})

target_link_libraries(gam-create ${ZLIB_LIBRARIES})
target_link_libraries(gam-create ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gam-create BamTools) #target_link_libraries(gam-ngs ${PROJECT_SOURCE_DIR}/lib/BamFile/libbamtools.a)
target_link_libraries(gam-create ${Boost_LIBRARIES})

# GAM-MERGE executable
add_executable(gam-merge src/gam-merge.cc src/Merge.cc src/Options.cc src/OptionsMerge.cc ${GAMNGSLIB_SRC_FILES})

target_link_libraries(gam-merge ${ZLIB_LIBRARIES})
target_link_libraries(gam-merge ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gam-merge BamTools) #target_link_libraries(gam-ngs ${PROJECT_SOURCE_DIR}/lib/BamFile/libbamtools.a)
target_link_libraries(gam-merge ${Boost_LIBRARIES})

# GAM-N50 executable
add_executable(gam-n50 src/n50.cc)
