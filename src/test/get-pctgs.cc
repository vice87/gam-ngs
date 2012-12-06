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

#include "Merge.hpp"

#include <iostream>

#include <map>
#include <list>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

#include "assembly/contig.hpp"
#include "assembly/io_contig.hpp"

#include "pool/HashContigMemPool.hpp"
#include "UtilityFunctions.hpp"

int main(int argc, char *argv[])
{
	if( argc != 3 )
	{
		std::cerr << "usage: " << argv[0] << "<GAM-assembly.fasta> <pctgs file>" << std::endl;
		return 1;
	}

	std::string fasta_file = argv[1];
	std::string pctgs_file = argv[2];

	std::map< std::string, uint64_t > pnames;

	std::string line;
	std::string pctg_name, trash;

	std::ifstream pfs( pctgs_file.c_str() );

	getline( pfs, line ); // discard first line

	do
	{
		getline( pfs, line ); // read filename

		std::stringstream ss(line);
		ss >> pctg_name >> trash >> trash >> trash >> trash >> trash;

		if( pnames.find( pctg_name ) == pnames.end() ) pnames[pctg_name] = 1; else pnames[pctg_name] += 1;
	}
	while( line != "" );

	pfs.close();

	// load only contigs
	HashContigMemPool pctg_pool;
	std::string outfile = fasta_file + std::string(".pctgs.fasta");

	std::ifstream ifs( fasta_file.c_str(), std::ifstream::in );
	std::ofstream ofs( outfile.c_str(), std::ios::out );

	while( !ifs.eof() )
	{
		std::string ctg_name;
		pctg_pool.readNextContigID( ifs, ctg_name );

		Contig ctg;
		pctg_pool.readNextSequence( ifs, ctg );

		if( pnames.find( ctg_name ) != pnames.end() )
		{
			if( pnames[ctg_name] > 1 ) ofs << ctg << std::endl;
		}
	}

	ifs.close();
	ofs.close();

	return 0;
}