/*
 *  This file is part of rNA.
 *  Copyright (c) 2011 by Cristian Del Fabbro <delfabbro@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Alexandru Tomescu <alexandru.tomescu@uniud.it>, and
 *  Alberto Policriti <policriti@uniud.it>
 *
 *   rNA is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rNA is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with rNA.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "CreateBlocks.hpp"
#include "OptionsCreate.hpp"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include <google/sparse_hash_map>
#include <boost/filesystem.hpp>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "bam/MultiBamReader.hpp"
#include "assembly/Read.hpp"
#include "assembly/Block.hpp"
#include "UtilityFunctions.hpp"

using namespace BamTools;
using google::sparse_hash_map;

extern OptionsCreate g_options;

namespace modules
{

void CreateBlocks::execute()
{
	struct stat st;

	time_t t1 = time(NULL);

	std::cout << "[main] opening BAM files" << std::endl;

	// load master BAM filenames and min/max insert sizes
	std::vector< std::string > masterBamFiles; // vector of master BAM filenames
	std::vector< int32_t > master_minInsert, master_maxInsert;
	loadBamFileNames( g_options.masterBamFile, masterBamFiles, master_minInsert, master_maxInsert );

	// check master's BAM files existence
	for( int i=0; i < masterBamFiles.size(); i++ )
	{
		boost::filesystem::path p(masterBamFiles[i].c_str());
		if( !boost::filesystem::exists(p) || !boost::filesystem::is_regular_file(p) )
		{
			std::cerr << "[error] master BAM file \"" << masterBamFiles[i] << "\" doesn't exist" << std::endl;
			exit(1);
		}
	}

	// load slaves BAM filenames and min/max insert sizes
	std::vector< std::string > slaveBamFiles;
	std::vector< int32_t > slave_minInsert, slave_maxInsert;

	loadBamFileNames( g_options.slaveBamFile, slaveBamFiles, slave_minInsert, slave_maxInsert );

	// check slaves' BAM files existence
	for( int i=0; i < slaveBamFiles.size(); i++ )
	{
		boost::filesystem::path p(slaveBamFiles[i].c_str());
		if( !boost::filesystem::exists(p) || !boost::filesystem::is_regular_file(p) )
		{
			std::cerr << "[error] slave BAM file \"" << slaveBamFiles[i] << "\" doesn't exist" << std::endl;
			exit(1);
		}
	}

	/* OPEN MASTER BAM AND LOAD READS IN MEMORY */

	MultiBamReader masterBam; // master (multi) BAM reader
	masterBam.Open( masterBamFiles ); // open master BAM files
	masterBam.setMinMaxInsertSizes( master_minInsert, master_maxInsert );

	std::cout << "[main] loading reads in memory" << std::endl;

	std::vector< std::vector<uint32_t> > masterCoverage;
	sparse_hash_map< std::string, Read > masterReadMap_1, masterReadMap_2;

	// load uniquely mapped reads of the master, while updating master contig's coverage and inserts stats
	Read::loadReadsMap( masterBam, masterReadMap_1, masterReadMap_2, masterCoverage );

	// output inserts statistics for master assembly
	std::string isize_stats_file = g_options.masterBamFile + ".isize";
	masterBam.writeStatsToFile( isize_stats_file );

	masterBam.Close(); // close master bam (no longer needed)

	time_t t2 = time(NULL);
	std::cout << "[main] reads loaded in " << formatTime(t2-t1) << std::endl;

	std::cout << "[main] finding blocks" << std::endl;

	std::vector<Block> blocks;

	/* OPEN SLAVE BAM AND BUILD BLOCKS */

	MultiBamReader slaveBam; // slave (multi) BAM reader
	slaveBam.Open( slaveBamFiles ); // open slave BAM files
	slaveBam.setMinMaxInsertSizes( slave_minInsert, slave_maxInsert );

	std::vector< std::vector<uint32_t> > slaveCoverage;

	// build blocks, compute slave contig's coverage and inserts stats
	Block::findBlocks( blocks, slaveBam, g_options.minBlockSize,
					   masterReadMap_1, masterReadMap_2, slaveCoverage );


	/* COMPUTE COVERAGE OF THE BLOCKS */
	Block::updateCoverages( blocks, masterCoverage, slaveCoverage );

	// output inserts statistcs for slave assembly
	isize_stats_file = g_options.slaveBamFile + ".isize";
	slaveBam.writeStatsToFile( isize_stats_file );

	slaveBam.Close(); // close current slave (no longer needed)

	std::cout << "[main] blocks found = " << blocks.size() << std::endl;

	std::cout << "[main] writing blocks on file: " << getPathBaseName( g_options.outputFilePrefix ) << std::endl;
	Block::writeBlocks( g_options.outputFilePrefix + ".blocks", blocks );

	std::cout << "[main] total execution time = " << formatTime( time(NULL)-t1 ) << std::endl;
}

} // namespace modules
