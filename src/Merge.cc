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

#include <boost/graph/graphviz.hpp>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "types.hpp"
#include "assembly/Read.hpp"
#include "assembly/Block.hpp"
#include "graphs/PairingEvidencesGraph.hpp"
#include "pctg/PairedContig.hpp"
#include "pool/HashContigMemPool.hpp"
#include "OrderingFunctions.hpp"
#include "PartitionFunctions.hpp"
#include "ThreadedBuildPctg.code.hpp"
#include "UtilityFunctions.hpp"

using namespace BamTools;

//extern pthread_mutex_t g_badAlignMutex;
//extern std::vector<uint64_t> ext_fpi;
//extern std::vector<uint64_t> ext_fpi_2;
//extern std::ofstream g_badAlignStream;
//extern std::ofstream ext_ba_desc_stream;
//extern uint64_t ext_ba_count;

namespace modules
{

void Merge::execute(const Options &options)
{
	//pthread_mutex_init(&g_badAlignMutex,NULL);
	//ext_fpi = std::vector<uint64_t>(1000,0);
	//ext_fpi_2 = std::vector<uint64_t>(1000,0);
	//g_badAlignStream.open( (options.outputFilePrefix + ".bad_align").c_str() );
	//ext_ba_desc_stream.open( (options.outputFilePrefix + ".ba_desc").c_str() );
	//ext_ba_count = 1;

	struct stat st;
	time_t tStart = time(NULL);

	BamReader bamReader;

	std::cout << "[loading blocks]" << std::endl;
	std::vector<Block> blocks = Block::readBlocks( options.blocksFile, options.minBlockSize );
	std::cout << "[blocks loaded: " << blocks.size() << "]" << std::endl;

	std::cout << "[loading reference data]" << std::endl;

	// load master bam sequences
	bamReader.Open( options.masterBamFile );
	BamTools::RefVector mcRef = bamReader.GetReferenceData();

	// load slaves bams sequences
	std::vector< BamTools::RefVector > scRef( options.slaveBamFiles.size() );

	for( int i=0; i < options.slaveBamFiles.size(); i++ )
	{
		bamReader.Open( options.slaveBamFiles[i] );
		scRef[i] = bamReader.GetReferenceData();
	}

	bamReader.Close();

	// keep only blocks between contigs that share at least minEvidence blocks.
	//std::vector< Block > outBlocks = filterBlocksByPairingEvidences( blocks, minEvidence );

	// remove adjacent blocks if their frames overlap
	// blocks = Block::filterBlocksByOverlaps( blocks );

	std::cout << "[filtering " << blocks.size() << " blocks by coverage]" << std::endl;
	blocks = Block::filterBlocksByCoverage( blocks, 0.75 );

	std::cout << "[filtering " << blocks.size() << " blocks by length]" << std::endl;
	blocks = Block::filterBlocksByLength( blocks, mcRef, scRef, 500 );

	std::cout << "[remaining blocks: " << blocks.size() << "]" << std::endl;

	// create the graph of assemblies, remove cycles and keep remaining blocks.
	std::cout << "[removing cycles from assemblies graph]" << std::endl;
	std::list< std::vector<Block> > pcblocks = partitionBlocks( blocks, options );

	std::cout << "[loading contigs lengths]" << std::endl;

	std::map< std::string, int32_t > masterContigs;
	std::vector< std::map<std::string,int32_t> > slaveCtgsVect( options.slaveBamFiles.size() );

	BamTools::RefVector::const_iterator ref_iter;

	for( ref_iter = mcRef.begin(); ref_iter != mcRef.end(); ref_iter++ )
		masterContigs[ ref_iter->RefName ] = ref_iter->RefLength;

	for( int i=0; i < slaveCtgsVect.size(); i++ )
	{
		for( ref_iter = scRef[i].begin(); ref_iter != scRef[i].end(); ref_iter++ )
			slaveCtgsVect[i][ ref_iter->RefName ] = ref_iter->RefLength;
	}

	std::cout << "[loading contigs in memory]" << std::endl;

	ExtContigMemPool masterPool, slavePool(scRef.size());
	HashContigMemPool pctgPool;
	// load master and slave contigs in memory
	masterPool.loadPool( 0, options.masterFastaFile, masterContigs );
	for( int i=0; i < slaveCtgsVect.size(); i++ ) slavePool.loadPool( i, options.slaveFastaFiles[i], slaveCtgsVect[i] );

	std::cout << "[building paired contigs]" << std::endl;
	// build paired contigs (threaded)
	ThreadedBuildPctg tbp( pcblocks, &pctgPool, &masterPool, &slavePool, &mcRef, &scRef, options );
	std::pair< std::list<PairedContig>, std::vector<bool> > result = tbp.run();

	IdType pctgNum( result.first.size() );
	std::cout << "[paired contigs built: " << pctgNum << "]" << std::endl;

	// output assemblies made exclusively by contigs involved in merging
	/*std::fstream masterMergeFile( (options.outputFilePrefix + ".onlymaster.fasta").c_str(), std::fstream::out );
	std::fstream slaveMergeFile( (options.outputFilePrefix + ".onlyslave.fasta").c_str(), std::fstream::out );
	std::fstream linearMergeFile( (options.outputFilePrefix + ".onlymerge.fasta").c_str(), std::fstream::out );

	std::list<PairedContig>::iterator pctg;
	for( pctg = result.first.begin(); pctg != result.first.end(); pctg++ )
	{
		linearMergeFile << *pctg << std::endl;

		std::map< std::pair<IdType,IdType>, ContigInPctgInfo >::const_iterator ctg_id;
		for( ctg_id = (pctg->getMasterCtgMap()).begin(); ctg_id != (pctg->getMasterCtgMap()).end(); ctg_id++ )
			masterMergeFile << masterPool.get( (ctg_id->first).first, mcRef[(ctg_id->first).second].RefName ) << std::endl;

		for( ctg_id = (pctg->getSlaveCtgMap()).begin(); ctg_id != (pctg->getSlaveCtgMap()).end(); ctg_id++ )
			slaveMergeFile << slavePool.get( (ctg_id->first).first, scRef[(ctg_id->first).first][(ctg_id->first).second].RefName ) << std::endl;
	}

	masterMergeFile.close();
	slaveMergeFile.close();
	linearMergeFile.close();*/

	slavePool.resize(0); // slaves pool is no longer needed

	std::list<IdType> ctgIds;
	for(IdType i = 0; i < result.second.size(); i++) if( !result.second.at(i) ) ctgIds.push_back(i);

	// add unused contigs in paired contigs pool
	generateSingleCtgPctgs( result.first, ctgIds, &masterPool, &mcRef, pctgNum);

	// save paired contig pool to disk
	std::cout << "[writing paired contigs on file: " << ( options.outputFilePrefix + ".fasta" ) << "]" << std::endl << std::flush;

	result.first.sort( orderPctgsByName );

	std::ofstream outFasta((options.outputFilePrefix + ".fasta").c_str(),std::ios::out);
	std::list< PairedContig >::const_iterator i;
	for( i = result.first.begin(); i != result.first.end(); i++ ) outFasta << *i << std::endl;
	outFasta.close();

	// save paired contigs descriptors to file
	std::cout << "[writing paired contigs descriptors on file: " << ( options.outputFilePrefix + ".pctgs" ) << "]" << std::endl << std::flush;
	std::fstream pctgDescFile( (options.outputFilePrefix + ".pctgs").c_str(), std::fstream::out );
	writePctgDescriptors( pctgDescFile, result.first, mcRef, scRef );
	pctgDescFile.close();

	// save IDs of (slave) contigs NOT merged
	std::cout << "[writing unused slave contigs on file: " << ( options.outputFilePrefix + ".unused" ) << "]" << std::endl << std::flush;
	std::fstream unusedCtgsFile( (options.outputFilePrefix + ".unused").c_str(), std::fstream::out );
	std::vector< std::vector<bool> > usedCtgs;
	for( size_t i=0; i < scRef.size(); i++ ) usedCtgs.push_back( std::vector<bool>( scRef[i].size(), false ) );

	for( std::list< PairedContig >::const_iterator pctg = result.first.begin(); pctg != result.first.end(); pctg++ )
	{
		typedef std::map< std::pair<IdType,IdType>, ContigInPctgInfo > ContigInfoMap;
		ContigInfoMap slaveCtgs = pctg->getSlaveCtgMap();

		for( ContigInfoMap::const_iterator ctg = slaveCtgs.begin(); ctg != slaveCtgs.end(); ctg++ )
			usedCtgs[(ctg->first).first][(ctg->first).second] = true;
	}

	for( unsigned int i = 0; i < usedCtgs.size(); i++ )
	{
		for( size_t j=0; j < usedCtgs[i].size(); j++ )
			if( !usedCtgs[i][j] ){ unusedCtgsFile << i << "\t" << scRef[i][j].RefName << "\n"; }
	}

	unusedCtgsFile.close();

	std::cout << "[total execution time = " << formatTime( time(NULL)-tStart ) << "]" << std::endl;

	/*std::ofstream ofs_fpi( (options.outputFilePrefix + ".fpi").c_str() );
	int max_frame_num = 0;
	for( int i=1; i < ext_fpi.size(); i++ ){ if( (ext_fpi[i] > 0 || ext_fpi_2[i] > 0) && i > max_frame_num) max_frame_num = i; }
	ofs_fpi << "          \t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << i << "\t"; } ofs_fpi << std::endl;
	ofs_fpi << "successful\t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << ext_fpi[i] << "\t"; } ofs_fpi << std::endl;
	ofs_fpi << "failed    \t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << ext_fpi_2[i] << "\t"; } ofs_fpi << std::endl;
	ofs_fpi.close();

	g_badAlignStream.close();*/
}

} // namespace modules
