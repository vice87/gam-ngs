/*
 * File:   hash_test.cc
 * Author: vice
 *
 * Created on March 26, 2012, 6:33 PM
 */

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>

#include <time.h>

#include "assembly/Frame.hpp"
#include "assembly/Block.hpp"
#include "assembly/Read.hpp"

int main(int argc, char** argv)
{
	if( argc != 2 )
	{
		std::cout << "Usage: blocks-coverage <file.blocks>" << std::endl;
		return 1;
	}

	std::string blocksFile = argv[1];
	std::vector<Block> blocks = Block::readBlocks( blocksFile );

	std::vector< uint64_t > coverage(51,0);

	uint32_t mcRatio, scRatio, maxRatio;

	for( size_t i=0; i < blocks.size(); i++ )
	{
		Frame& mf = blocks[i].getMasterFrame();
		Frame& sf = blocks[i].getSlaveFrame();

		mcRatio = (mf.getReadsLen() == 0) ? 0 : (100 * mf.getBlockReadsLen()) / mf.getReadsLen();
		scRatio = (sf.getReadsLen() == 0) ? 0 : (100 * sf.getBlockReadsLen()) / sf.getReadsLen();
		maxRatio = std::max(mcRatio,scRatio);

		if( maxRatio > 100 )
		{
			std::cerr << "copertura maggiore del 100%, c'è un errore!" << std::endl;
			continue;
		}

		uint32_t idx = maxRatio/2;
		coverage[idx]++;
	}

	std::string outputFile = blocksFile + ".hist";

	std::ofstream ofs( outputFile.c_str() );
	for( size_t i=0; i < coverage.size(); i++ ) ofs << i*2 << "\t" << coverage[i] << "\n";
	ofs.close();

    return 0;
}

