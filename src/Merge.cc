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
#include "OptionsMerge.hpp"

#include "assembly/Block.hpp"
#include "assembly/Read.hpp"
#include "assembly/RefSequence.hpp"
#include "assembly/io_contig.hpp"
#include "bam/MultiBamReader.hpp"
#include "graphs/PairingEvidencesGraph.hpp"
#include "pctg/PairedContig.hpp"
#include "pctg/ThreadedBuildPctg.hpp"
#include "pctg/BuildPctgFunctions.hpp"
#include "OrderingFunctions.hpp"
#include "PartitionFunctions.hpp"
#include "UtilityFunctions.hpp"

using namespace BamTools;

extern OptionsMerge g_options;

extern std::ofstream g_badAlignStream;
//extern std::ofstream ext_ba_desc_stream;
//extern uint64_t ext_ba_count;

extern std::vector<uint64_t> ext_fpi;
extern std::vector<uint64_t> ext_fpi_2;
extern pthread_mutex_t g_badAlignMutex;

std::ofstream _g_statsFile;

MultiBamReader masterBam;
MultiBamReader masterMpBam;
MultiBamReader slaveBam;
MultiBamReader slaveMpBam;

namespace modules
{

void Merge::execute()
{
	pthread_mutex_init(&g_badAlignMutex,NULL);
	ext_fpi = std::vector<uint64_t>(1000,0);
	ext_fpi_2 = std::vector<uint64_t>(1000,0);

	g_badAlignStream.open( (g_options.outputFilePrefix + ".bad_align").c_str() );
	//ext_ba_desc_stream.open( (options.outputFilePrefix + ".ba_desc").c_str() );
	//ext_ba_count = 1;

	struct stat st;
	time_t tStart = time(NULL);

	_g_statsFile.open( (g_options.outputFilePrefix + ".stats").c_str(), std::ios::out ); // open statistics (output) file

	std::list<Block> blocks;

	std::cout << "[main] loading blocks" << std::endl;
	Block::loadBlocks( g_options.blocksFile, blocks, g_options.minBlockSize );
	std::cout << "[main] loaded blocks = " << blocks.size() << std::endl;

	std::cout << "[main] loading BAMs data" << std::endl;

	std::vector< int32_t > minInsert, maxInsert;

	/* OPEN MASTER BAM FILES */

	std::vector< std::string > masterBamFiles; // vector of master BAM filenames
	loadBamFileNames( g_options.masterBamFile, masterBamFiles, minInsert, maxInsert ); // load master BAM alignments filenames

	masterBam.Open( masterBamFiles ); // open master BAM files
	if( masterBam.size() == 0 )
	{
		std::cerr << "[bam] error: cannot open any master's PE-alignment specified in " << g_options.masterBamFile << std::endl;
		exit(1);
	}
	
	masterBam.setMinMaxInsertSizes( minInsert, maxInsert );

	// open inserts statistics for master alignments
	if( g_options.masterISizeFile != "" )
	{
		if( stat(g_options.masterISizeFile.c_str(),&st) != 0 ) // if statistics file do not exist, create it
		{
			std::cout << "[bam] computing statistics of master's PE-alignments" << std::endl;
			masterBam.computeStatistics();
			masterBam.writeStatsToFile( g_options.masterISizeFile );
		}
		
		masterBam.readStatsFromFile( g_options.masterISizeFile ); 
	}

	std::cout << "[bam] master PE-alignments file " << getPathBaseName(g_options.masterBamFile) << " successfully opened:" << std::endl;
	for( size_t i=0; i < masterBam.size(); i++ ) 
		std::cout << "      " << masterBam[i].GetFilename() 
			<< "\n         inserts size = " << masterBam.getISizeMean(i) << " +/- " << masterBam.getISizeStd(i) 
			<< "\tcoverage = " << masterBam.getCoverage(i) << std::endl;

	/* OPEN MASTER MP BAM FILES */

	if( g_options.masterMpBamFile != "" )
	{
		std::vector< std::string > masterMpBamFiles; // vector of master BAM filenames
		loadBamFileNames( g_options.masterMpBamFile, masterMpBamFiles, minInsert, maxInsert ); // load master BAM alignments filenames

		masterMpBam.Open( masterMpBamFiles ); // open master BAM files
		if( masterMpBam.size() == 0 )
		{
			std::cerr << "[bam] error: cannot open any master's MP-alignment specified in " << g_options.masterMpBamFile << std::endl;
			exit(1);
		}
		
		masterMpBam.setMinMaxInsertSizes( minInsert, maxInsert );

		if( g_options.masterMpISizeFile != "" ) 
		{
			if( stat(g_options.masterMpISizeFile.c_str(),&st) != 0 ) // if statistics file do not exist, create it
			{
				std::cout << "[bam] computing statistics of master's MP-alignments" << std::endl;
				masterMpBam.computeStatistics();
				masterMpBam.writeStatsToFile( g_options.masterMpISizeFile );
			}
			
			masterMpBam.readStatsFromFile( g_options.masterMpISizeFile ); // open inserts statistics
		}

		std::cout << "[bam] master MP-alignments file " << getPathBaseName(g_options.masterMpBamFile) << " successfully opened:" << std::endl;
		for( size_t i=0; i < masterMpBam.size(); i++ ) std::cout << "      " << masterMpBam[i].GetFilename() 
			<< "\n         inserts size = " << masterMpBam.getISizeMean(i) << " +/- " << masterMpBam.getISizeStd(i) 
			<< "\tcoverage = " << masterMpBam.getCoverage(i) << std::endl;
	}

	/* OPEN SLAVE BAM FILES */

	std::vector< std::string > slaveBamFiles; // vector of slave BAM filenames
	loadBamFileNames( g_options.slaveBamFile, slaveBamFiles, minInsert, maxInsert ); // load slaves BAM alignments filenames

	slaveBam.Open( slaveBamFiles );
	if( slaveBam.size() == 0 )
	{
		std::cerr << "[bam] error: cannot open any slave's PE-alignment specified in " << g_options.slaveBamFile << std::endl;
		exit(1);
	}
	
	slaveBam.setMinMaxInsertSizes( minInsert, maxInsert );

	if( g_options.slaveISizeFile != "" ) 
	{
		if( stat(g_options.slaveISizeFile.c_str(),&st) != 0 ) // if statistics file do not exist, create it
		{
			std::cout << "[bam] computing statistics of slave's PE-alignments" << std::endl;
			slaveBam.computeStatistics();
			slaveBam.writeStatsToFile( g_options.slaveISizeFile );
		}
		
		slaveBam.readStatsFromFile( g_options.slaveISizeFile ); // open inserts statistics
	}

	std::cout << "[bam] slave PE-alignments file " << getPathBaseName(g_options.slaveBamFile) << " successfully opened:" << std::endl;
	for( size_t i=0; i < slaveBam.size(); i++ ) std::cout << "      "<< slaveBam[i].GetFilename() 
		<< "\n         inserts size = " << slaveBam.getISizeMean(i) << " +/- " << slaveBam.getISizeStd(i)
		<< "\tcoverage = " << slaveBam.getCoverage(i) << std::endl;

	/* OPEN SLAVE MP BAM FILES */

	if( g_options.slaveMpBamFile != "" )
	{
		std::vector< std::string > slaveMpBamFiles; // vector of slave BAM filenames
		loadBamFileNames( g_options.slaveMpBamFile, slaveMpBamFiles, minInsert, maxInsert ); // load slaves BAM alignments filenames

		slaveMpBam.Open( slaveMpBamFiles );
		if( slaveMpBam.size() == 0 )
		{
			std::cerr << "[bam] error: cannot open any slave's MP-alignment specified in " << g_options.slaveMpBamFile << std::endl;
			exit(1);
		}
		
		slaveMpBam.setMinMaxInsertSizes( minInsert, maxInsert );

		if( g_options.slaveMpISizeFile != "" ) 
		{
			if( stat(g_options.slaveMpISizeFile.c_str(),&st) != 0 ) // if statistics file do not exist, create it
			{
				std::cout << "[bam] computing statistics of slave's MP-alignments" << std::endl;
				slaveMpBam.computeStatistics();
				slaveMpBam.writeStatsToFile( g_options.slaveMpISizeFile );
			}
			
			slaveMpBam.readStatsFromFile( g_options.slaveMpISizeFile ); // open inserts statistics
		}

		std::cout << "[bam] slave MP-alignments file " << getPathBaseName(g_options.slaveMpBamFile) << " successfully opened:" << std::endl;
		for( size_t i=0; i < slaveMpBam.size(); i++ ) std::cout << "      "<< slaveMpBam[i].GetFilename() 
			<< "\n         inserts size = " << slaveMpBam.getISizeMean(i) << " +/- " << slaveMpBam.getISizeStd(i)
			<< "\tcoverage = " << slaveMpBam.getCoverage(i) << std::endl;
	}

	/* LOAD SEQUENCES DATA */

	std::cout << "[main] loading contigs data" << std::endl;

	uint64_t master_ctgs = masterBam.GetReferenceData().size();
	uint64_t slave_ctgs = slaveBam.GetReferenceData().size();

	RefSequence masterRef(master_ctgs);
	for( uint64_t i=0; i < master_ctgs; i++ )
	{
		masterRef[i].RefName = masterBam.GetReferenceData().at(i).RefName;
		masterRef[i].RefLength = masterBam.GetReferenceData().at(i).RefLength;
	}

	RefSequence slaveRef(slave_ctgs);
	for( uint64_t i=0; i < slave_ctgs; i++ )
	{
		slaveRef[i].RefName = slaveBam.GetReferenceData().at(i).RefName;
		slaveRef[i].RefLength = slaveBam.GetReferenceData().at(i).RefLength;
	}

	/* BLOCKS FILTERING */

	std::set< std::pair<int32_t,int32_t> > sl_blocks;
	getSingleLinkBlocks( blocks, sl_blocks );

	std::set<int32_t> masterNBC_BF, slaveNBC_BF; // contigs without blocks (before filtering)
	getNoBlocksContigs( masterRef, slaveRef, blocks, masterNBC_BF, slaveNBC_BF );
	
	double min_cov = std::min( masterBam.getGlobCoverage(), slaveBam.getGlobCoverage() ) / 2;

	std::cout << "[main] filtering blocks by coverage" << std::endl;
	Block::filterBlocksByCoverage( blocks, sl_blocks, min_cov, g_options.coverageThreshold );
	std::cout << "[main] remaining blocks = " << blocks.size() << std::endl;
	
	//TODO: sistemare la funzione per eliminare blocchi dovuti a repeats
	//Block::filterBlocksByOverlaps(blocks);

	//TODO: migliorare il filtraggio dei blocchi per lunghezza
	//std::cout << "[main] filtering blocks by length" << std::endl;
	//Block::filterBlocksByLength( blocks, masterRef, slaveRef, sl_blocks, 500 );
	//std::cout << "[main] length filtered blocks = " << blocks.size() << std::endl;

	std::set<int32_t> masterNBC_AF, slaveNBC_AF; // contigs without blocks (after filtering)
	getNoBlocksAfterFilterContigs( masterRef, slaveRef, blocks, masterNBC_BF, slaveNBC_BF, masterNBC_AF, slaveNBC_AF );

	/* PARTITIONING BLOCKS */

	std::cout << "[main] partitioning blocks" << std::endl;
	std::list< std::vector<Block> > pcblocks = partitionBlocks( blocks );

	/* LOADING CONTIGS SEQUENCES IN MEMORY */

	std::cout << "[main] loading contigs sequences" << std::endl;

	std::map< std::string, int32_t > masterCtg2Id; // master contig Name => ID
	std::map< std::string, int32_t > slaveCtg2Id;  // slave contig Name => ID

	for( size_t i=0; i < masterRef.size(); i++ ) masterCtg2Id[ masterRef[i].RefName ] = i;
	for( size_t i=0; i < slaveRef.size(); i++ ) slaveCtg2Id[ slaveRef[i].RefName ] = i;

	loadSequences( g_options.masterFastaFile, masterRef, masterCtg2Id );
	loadSequences( g_options.slaveFastaFile, slaveRef, slaveCtg2Id );

	// output slave contigs with no blocks (before filtering)
	std::string noblocks_fasta_file = g_options.outputFilePrefix + ".noblocks.BF.fasta";
	std::cout << "[merge] writing contigs with no blocks to file \"" << noblocks_fasta_file << "\"" << std::endl;
	std::ofstream noblocks_fasta_stream( noblocks_fasta_file.c_str() );
	for( std::set< int32_t >::const_iterator it = slaveNBC_BF.begin(); it != slaveNBC_BF.end(); it++ )
	{
		const Contig *ctg = slaveRef[*it].Sequence;
		noblocks_fasta_stream << *ctg << std::endl;
	}
	noblocks_fasta_stream.close();
	
	// output slave contigs with no blocks (after filtering)
	noblocks_fasta_file = g_options.outputFilePrefix + ".noblocks.AF.fasta";
	std::cout << "[merge] writing contigs with no blocks (after filtering) to file \"" << noblocks_fasta_file << "\"" << std::endl;
	noblocks_fasta_stream.open( noblocks_fasta_file.c_str() );
	for( std::set< int32_t >::const_iterator it = slaveNBC_AF.begin(); it != slaveNBC_AF.end(); it++ )
	{
		const Contig *ctg = slaveRef[*it].Sequence;
		noblocks_fasta_stream << *ctg << std::endl;
	}
	noblocks_fasta_stream.close();
	
	// build paired contigs (threaded)
	std::cout << "[merge] building paired contigs" << std::endl;
	ThreadedBuildPctg tbp( pcblocks, masterRef, slaveRef );
	std::list<PairedContig> *result = tbp.run();
	
	// assign unique IDs to paired contigs
	uint64_t pctg_id = 0;
	for( std::list<PairedContig>::iterator pctg = result->begin(); pctg != result->end(); pctg++ )
	{
		pctg->setId(pctg_id);
		pctg_id++;
	}
	
	std::cout << "[merge] paired contigs built = " << pctg_id << std::endl;
	
	// TODO: sistemare codice commentato qui sotto
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

	// TODO: sistemare codice commentato qui sotto
	// save IDs of (slave) contigs NOT merged
	/*std::cout << "[merge] writing not merged (slave) contigs on file \"" << ( g_options.outputFilePrefix + ".notmerged.fasta" ) << "\"" << std::endl;
	std::fstream unusedCtgsFile( (g_options.outputFilePrefix + ".notmerged.fasta").c_str(), std::fstream::out );
	std::vector<bool> usedCtgs( slaveRef.size(), false );

	for( std::list< PairedContig >::const_iterator pctg = result->begin(); pctg != result->end(); pctg++ )
	{
		typedef std::map< int32_t, ContigInPctgInfo > ContigInfoMap;
		const ContigInfoMap& slaveCtgs = pctg->getSlaveCtgMap();

		for( ContigInfoMap::const_iterator ctg = slaveCtgs.begin(); ctg != slaveCtgs.end(); ctg++ )
			usedCtgs[ctg->first] = true;
	}

	for( std::set<int32_t>::const_iterator it = slaveNBC_BF.begin(); it != slaveNBC_BF.end(); it++ )
		usedCtgs[*it] = true;

	for( std::set<int32_t>::const_iterator it = slaveNBC_AF.begin(); it != slaveNBC_AF.end(); it++ )
		usedCtgs[*it] = true;

	for( size_t i=0; i < usedCtgs.size(); i++ )
		if( !usedCtgs[i] ) unusedCtgsFile << *(slaveRef[i].Sequence) << "\n";

	unusedCtgsFile.close();
	usedCtgs.clear(); // linea commentata perchè mi servono dopo per le regioni duplicate
	*/

	// delete slave contigs, which are no longer needed
	for( size_t i=0; i < slaveRef.size(); i++ ) delete slaveRef[i].Sequence;
	slaveNBC_BF.clear(); slaveNBC_AF.clear();

	// get master contigs that have not been included in any pctg
	std::vector< bool > usedMasterCtgs( masterRef.size(), false );
	for( std::list< PairedContig >::iterator pctg = result->begin(); pctg != result->end(); pctg++ )
	{
		std::set<int32_t>::const_iterator it;
		for( it = (pctg->getMasterCtgIdSet()).begin(); it != (pctg->getMasterCtgIdSet()).end(); it++ ) usedMasterCtgs[*it] = true;
	}

	std::list<int32_t> ctgIds;
	for(int32_t i=0; i < usedMasterCtgs.size(); i++) if( !usedMasterCtgs[i] ) ctgIds.push_back(i);

	// add unused contigs in paired contigs pool
	uint64_t old_pctg_id = pctg_id;
	generateSingleCtgPctgs( *result, ctgIds, masterRef, pctg_id);

	// save paired contig pool to disk
	std::cout << "[merge] writing paired contigs on file \"" << ( g_options.outputFilePrefix + ".gam.fasta" ) << "\"" << std::endl;

	std::ofstream outFasta((g_options.outputFilePrefix + ".gam.fasta").c_str(),std::ios::out);
	for( std::list< PairedContig >::const_iterator pctg = result->begin(); pctg != result->end(); pctg++ ) outFasta << *pctg << std::endl;
	outFasta.close();

	// save paired contigs descriptors to file
	std::cout << "[merge] writing paired contigs descriptors on file \"" << ( g_options.outputFilePrefix + ".pctgs" ) << "\"" << std::endl;
	std::fstream pctgDescFile( (g_options.outputFilePrefix + ".pctgs").c_str(), std::fstream::out );
	writePctgDescriptors( pctgDescFile, *result, masterRef, slaveRef, old_pctg_id );
	pctgDescFile.close();

	_g_statsFile.close();

	// DEBUG // TODO: incorporare meglio nel codice il calcolo delle seguenti statistiche
	/*std::ofstream ofs_fpi( (g_options.outputFilePrefix + ".fpi").c_str() );
	int max_frame_num = 0;
	for( int i=1; i < ext_fpi.size(); i++ ){ if( (ext_fpi[i] > 0 || ext_fpi_2[i] > 0) && i > max_frame_num) max_frame_num = i; }
	ofs_fpi << "          \t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << i << "\t"; } ofs_fpi << std::endl;
	ofs_fpi << "successful\t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << ext_fpi[i] << "\t"; } ofs_fpi << std::endl;
	ofs_fpi << "failed    \t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << ext_fpi_2[i] << "\t"; } ofs_fpi << std::endl;
	ofs_fpi.close();*/

	std::cout << "[merge] total execution time = " << formatTime( time(NULL)-tStart ) << std::endl;

	g_badAlignStream.close();
}

} // namespace modules
