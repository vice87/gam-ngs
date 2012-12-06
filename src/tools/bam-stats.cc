/*
 * File:   hash_test.cc
 * Author: vice
 *
 * Created on March 26, 2012, 6:33 PM
 */

#include <cstdlib>
#include <iostream>
#include <time.h>

#include "UtilityFunctions.hpp"

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "bam/MultiBamReader.hpp"

using namespace BamTools;

int main(int argc, char** argv)
{
	if( argc != 2 )
	{
		std::cerr << "Usage: " << argv[0] << " <BAM_list.txt>" << std::endl;
		return 1;
	}

	time_t start = time(NULL);

	std::string filename = argv[1];
	std::vector< std::string > bamFiles; // vector of BAM filenames
	std::vector< int32_t > minInsert, maxInsert; // min/max insert sizes to compute lib statistics
	MultiBamReader bamReader;

	loadBamFileNames( filename, bamFiles, minInsert, maxInsert );
	bamReader.Open( bamFiles ); // open master BAM files
	bamReader.setMinMaxInsertSizes( minInsert, maxInsert );

	BamAlignment align;
	while( bamReader.GetNextAlignment(align,true) );

	std::string isize_stats_file = filename + ".isize";
	bamReader.writeStatsToFile( isize_stats_file );

	bamReader.Close();

	time_t end = time(NULL);

	std::cout << "time = " << formatTime( end-start ) << std::endl;

    return 0;
}

