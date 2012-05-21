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

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <google/sparse_hash_map>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "bam/MultiBamReader.hpp"
#include "assembly/Read.hpp"
#include "assembly/Block.hpp"
#include "UtilityFunctions.hpp"

using namespace BamTools;
using google::sparse_hash_map;

namespace modules
{

void CreateBlocks::execute(const Options &options)
{
	struct stat st;
	time_t tStart = time(NULL);

	MultiBamReader masterBam; // master (multi) BAM reader

	std::vector< std::string > masterBamFiles; // vector of master BAM filenames
	loadFileNames( options.masterBamFiles, masterBamFiles );
	masterBam.Open( masterBamFiles ); // open master BAM files

	std::cout << "[loading reads in memory]" << std::endl;

	std::vector< std::vector<uint32_t> > masterCoverage;
	sparse_hash_map< std::string, Read > masterReadMap_1, masterReadMap_2;
	//masterReadMap.set_deleted_key("");

	// load only useful reads of the slave assembly
	Read::loadReadsMap( masterBam, masterReadMap_1, masterReadMap_2, masterCoverage );

	// output insert size statistics for master assembly
	std::string isize_stats_file = options.masterBamFiles + ".isize";
	masterBam.writeStatsToFile( isize_stats_file );

	time_t tEnd = time(NULL);
	std::cout << "[reads loaded in " << formatTime(tEnd-tStart) << "]" << std::endl;

	std::cout << "[finding blocks]" << std::endl;

	int32_t sid = 0;
	std::vector<Block> blocks;

	for( unsigned int i=0; i < options.slaveBamFiles.size(); i++ )
	{
		std::vector< std::string > slaveBamFiles; // vector of slave BAM filenames (of the current fp assembly)

		MultiBamReader slaveBam; // slave (multi) BAM reader
		loadFileNames( options.slaveBamFiles[i], slaveBamFiles );
		slaveBam.Open( slaveBamFiles ); // open current fp slave BAM files

		std::vector< std::vector<uint32_t> > slaveCoverage;
		Block::findBlocks( blocks, slaveBam, options.minBlockSize, masterReadMap_1, masterReadMap_2, slaveCoverage, sid );
		Block::updateCoverages( blocks, masterCoverage, slaveCoverage );

		// output insert size statistcs for each slave assembly
		isize_stats_file = options.masterBamFiles + ".isize";
		slaveBam.writeStatsToFile( isize_stats_file );

		sid++;
	}

	time_t tEnd2 = time(NULL);
	std::cout << "[" << blocks.size() << " blocks found in " << formatTime(tEnd2-tEnd) << "]" << std::endl;

	std::cout << "[writing blocks on file: " << getPathBaseName( options.outputFilePrefix ) << "]" << std::endl;
	Block::writeBlocks( options.outputFilePrefix + ".blocks", blocks );

	std::cout << "[total execution time = " << formatTime( time(NULL)-tStart ) << "]" << std::endl;
}

} // namespace modules
