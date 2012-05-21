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

#include "bam/MultiBamReader.hpp"

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

//extern std::ofstream g_badAlignStream;
//extern std::ofstream ext_ba_desc_stream;
//extern uint64_t ext_ba_count;

extern std::vector<uint64_t> ext_fpi;
extern std::vector<uint64_t> ext_fpi_2;
extern pthread_mutex_t g_badAlignMutex;
std::ofstream _g_statsFile;

namespace modules
{

void Merge::execute(const Options &options)
{
	pthread_mutex_init(&g_badAlignMutex,NULL);
	ext_fpi = std::vector<uint64_t>(1000,0);
	ext_fpi_2 = std::vector<uint64_t>(1000,0);

	//g_badAlignStream.open( (options.outputFilePrefix + ".bad_align").c_str() );
	//ext_ba_desc_stream.open( (options.outputFilePrefix + ".ba_desc").c_str() );
	//ext_ba_count = 1;

	struct stat st;
	time_t tStart = time(NULL);

	_g_statsFile.open( (options.outputFilePrefix + ".stats").c_str(), std::ios::out ); // open statistics (output) file

	// lunghezza media e dev.std. degli inserti dei BAM utilizzati
	// std::vector< double > masterInsertMean, masterInsertStd;
	// std::vector< std::vector<double> > slavesInsertMean, slavesInsertStd;

	std::cout << "[loading blocks]" << std::endl;
	std::vector<Block> blocks = Block::readBlocks( options.blocksFile, options.minBlockSize );
	std::cout << "[loaded blocks: " << blocks.size() << "]" << std::endl;

	std::cout << "[loading BAM data]\n" << std::endl;

	MultiBamReader masterBam; // master (multi) BAM reader
	std::vector< std::string > masterBamFiles; // vector of master BAM filenames
	loadFileNames( options.masterBamFiles, masterBamFiles ); // load master BAM alignments filenames

	masterBam.Open( masterBamFiles ); // open master BAM files
	if( options.masterISizeFiles != "" ) masterBam.readStatsFromFile( options.masterISizeFiles ); // open inserts statistics
	BamTools::RefVector mcRef = masterBam.GetReferenceData(); // load master bam sequences

	std::cout << options.masterBamFiles << ":" << std::endl;
	for( size_t i=0; i < masterBam.size(); i++ )
	{
		std::cout << masterBam[i].GetFilename() << std::endl;
		std::cout << "Insert stats => Mean = " << masterBam.getISizeMean(i) << "\tStd = " << masterBam.getISizeStd(i) << std::endl;
	}
	std::cout << std::endl;

	std::vector< MultiBamReader > slaveBams( options.slaveBamFiles.size() );
	std::vector< BamTools::RefVector > scRef( options.slaveBamFiles.size() ); // load slaves bams sequences

	for( int i=0; i < options.slaveBamFiles.size(); i++ )
	{
		std::vector< std::string > slaveBamFiles; // vector of slave BAM filenames (of the current fp assembly)
		loadFileNames( options.slaveBamFiles[i], slaveBamFiles ); // load slaves BAM alignments filenames

		slaveBams[i].Open( slaveBamFiles );
		if( options.slaveISizeFiles.size() > 0 ) slaveBams[i].readStatsFromFile( options.slaveISizeFiles[i] ); // open inserts statistics
		scRef[i] = slaveBams[i].GetReferenceData(); // load current slave bam sequences

		std::cout << options.slaveBamFiles[i] << ":" << std::endl;
		for( size_t j=0; j < slaveBams[i].size(); j++ )
		{
			std::cout << slaveBams[i][j].GetFilename() << std::endl;
			std::cout << "Insert stats => Mean = " << slaveBams[i].getISizeMean(j) << "\tStd = " << slaveBams[i].getISizeStd(j) << std::endl;
		}
		std::cout << std::endl;
	}

	// keep only blocks between contigs that share at least minEvidence blocks.
	//std::vector< Block > outBlocks = filterBlocksByPairingEvidences( blocks, minEvidence );

	// remove adjacent blocks if their frames overlap
	// blocks = Block::filterBlocksByOverlaps( blocks );

	std::cout << "[filtering blocks by coverage]" << std::endl;
	blocks = Block::filterBlocksByCoverage( blocks, 0.75 );
	std::cout << "[remaining blocks: " << blocks.size() << "]" << std::endl;

	std::cout << "[filtering blocks by length]" << std::endl;
	blocks = Block::filterBlocksByLength( blocks, mcRef, scRef, 500 );
	std::cout << "[remaining blocks: " << blocks.size() << "]" << std::endl;

	// create the graph of assemblies, remove cycles and keep remaining blocks.
	std::cout << "[removing cycles from assemblies graph]" << std::endl;
	std::list< std::vector<Block> > pcblocks = partitionBlocks( blocks, options );

	std::cout << "[loading contigs lengths]" << std::endl;

	std::map< std::string, int32_t > masterContigs;
	std::vector< std::map<std::string,int32_t> > slaveCtgsVect( options.slaveBamFiles.size() );

	BamTools::RefVector::const_iterator ref_iter;

	for( ref_iter = masterBam.GetReferenceData().begin(); ref_iter != masterBam.GetReferenceData().end(); ref_iter++ )
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
	ThreadedBuildPctg tbp( pcblocks, &pctgPool, &masterPool, &slavePool, masterBam, slaveBams, &mcRef, &scRef, options );
	std::pair< std::list<PairedContig>, std::vector<bool> > result = tbp.run();

	IdType pctgNum( result.first.size() );
	std::cout << "[paired contigs built: " << pctgNum << "]" << std::endl;


	// compute some statistics on merging and blocks' construction
	{
		std::set< int32_t > m_ctg_ids;
		std::vector< std::set< int32_t > > s_ctg_ids( options.slaveBamFiles.size() );

		uint64_t m_sum_len = 0, s_sum_len = 0;
		uint64_t m_ctg_num = 0, s_ctg_num = 0;

		int32_t aId, ctgId;

		uint64_t mb_ctgs = 0, sb_ctgs = 0;

		for( std::vector<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
		{
			const Frame& mf = b->getMasterFrame();
			const Frame& sf = b->getSlaveFrame();

			ctgId = mf.getContigId();

			if( m_ctg_ids.insert( ctgId ).second )
			{
				m_sum_len += mcRef[ctgId].RefLength;
				m_ctg_num++;
				mb_ctgs++;
			}

			aId = sf.getAssemblyId();
			ctgId = sf.getContigId();

			if( s_ctg_ids[aId].insert( ctgId ).second )
			{
				s_sum_len += scRef[aId][ctgId].RefLength;
				s_ctg_num++;
				sb_ctgs++;
			}
		}

		m_ctg_ids.clear();
		s_ctg_ids.clear();

		_g_statsFile << "[contigs without blocks]\n"
			<< "Master = " << masterPool.size() - mb_ctgs << " / " << masterPool.size() << "\n"
			<< "Slave  = " << slavePool.size() - sb_ctgs << " / " << slavePool.size() << "\n"
			<< "Master Total Length = " << m_sum_len << "\t Average = " << m_sum_len / double(m_ctg_num) << "\n"
			<< " Slave Total Length = " << s_sum_len << "\t Average = " << s_sum_len / double(s_ctg_num) << "\n"
			<< std::endl;

		std::set< int32_t > m_ctg_ids_1k;
		std::vector< std::set<int32_t> > s_ctg_ids_1k( options.slaveBamFiles.size() );

		uint64_t mb_ctgs_1k = 0, sb_ctgs_1k = 0, mp_1k = 0, sp_1k = 0;

		m_sum_len = s_sum_len = m_ctg_num = s_ctg_num = 0;

		for( std::vector<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
		{
			const Frame& mf = b->getMasterFrame();
			const Frame& sf = b->getSlaveFrame();

			ctgId = mf.getContigId();

			if( mcRef[ctgId].RefLength < 1000 ) continue;
			if( m_ctg_ids_1k.insert( ctgId ).second )
			{
				m_sum_len += mcRef.at(ctgId).RefLength;
				m_ctg_num++;
				mb_ctgs_1k++;
			}

			aId = sf.getAssemblyId();
			ctgId = sf.getContigId();

			if( scRef[aId][ctgId].RefLength < 1000 ) continue;
			if( s_ctg_ids_1k[aId].insert( ctgId ).second )
			{
				s_sum_len += scRef[aId][ctgId].RefLength;
				s_ctg_num++;
				sb_ctgs_1k++;
			}
		}

		m_ctg_ids_1k.clear();
		s_ctg_ids_1k.clear();

		for( size_t i=0; i < mcRef.size(); i++ ) if( mcRef[i].RefLength >= 1000 ) mp_1k++;
		for( size_t i=0; i < scRef.size(); i++ ) for( size_t j=0; j < scRef[i].size(); j++ ) if( scRef[i][j].RefLength >= 1000 ) sp_1k++;

		_g_statsFile << "[contigs >= 1K without blocks]\n"
		<< "Master = " << mp_1k - mb_ctgs_1k << " / " << mp_1k << "\n"
		<< "Slave  = " << sp_1k - sb_ctgs_1k << " / " << sp_1k << "\n"
		<< "Master Total Length = " << m_sum_len << "\t Average = " << m_sum_len / double(m_ctg_num) << "\n"
		<< " Slave Total Length = " << s_sum_len << "\t Average = " << s_sum_len / double(s_ctg_num) << "\n"
		<< std::endl;

		std::set< int32_t > m_ctg_ids_5k;
		std::vector< std::set<int32_t> > s_ctg_ids_5k( options.slaveBamFiles.size() );

		uint64_t mb_ctgs_5k = 0, sb_ctgs_5k = 0, mp_5k = 0, sp_5k = 0;
		m_sum_len = s_sum_len = m_ctg_num = s_ctg_num = 0;

		for( std::vector<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
		{
			const Frame& mf = b->getMasterFrame();
			const Frame& sf = b->getSlaveFrame();

			ctgId = mf.getContigId();

			if( mcRef[ctgId].RefLength < 5000 ) continue;
			if( m_ctg_ids_5k.insert( ctgId ).second )
			{
				m_sum_len += mcRef.at(ctgId).RefLength;
				m_ctg_num++;
				mb_ctgs_5k++;
			}

			aId = sf.getAssemblyId();
			ctgId = sf.getContigId();

			if( scRef[aId][ctgId].RefLength < 5000 ) continue;
			if( s_ctg_ids_5k[aId].insert( ctgId ).second )
			{
				s_sum_len += scRef[aId][ctgId].RefLength;
				s_ctg_num++;
				sb_ctgs_5k++;
			}
		}

		m_ctg_ids_5k.clear();
		s_ctg_ids_5k.clear();

		for( size_t i=0; i < mcRef.size(); i++ ) if( mcRef[i].RefLength >= 5000 ) mp_5k++;
		for( size_t i=0; i < scRef.size(); i++ ) for( size_t j=0; j < scRef[i].size(); j++ ) if( scRef[i][j].RefLength >= 5000 ) sp_5k++;

		_g_statsFile << "[contigs >= 5K without blocks]\n"
		<< "Master = " << mp_5k - mb_ctgs_5k << " / " << mp_5k << "\n"
		<< "Slave  = " << sp_5k - sb_ctgs_5k << " / " << sp_5k << "\n"
		<< "Master Total Length = " << m_sum_len << "\t Average = " << m_sum_len / double(m_ctg_num) << "\n"
		<< " Slave Total Length = " << s_sum_len << "\t Average = " << s_sum_len / double(s_ctg_num) << "\n"
		<< std::endl;

		_g_statsFile << "[merging stats]\n";

		uint64_t m_merged = 0, s_merged = 0, m_merged_tot = 0, s_merged_tot = 0;
		uint64_t merged_len = 0, merged_num = 0;
		std::list< uint32_t > m_merged_list, s_merged_list;

		// Per ogni PairedContig stampa il numero di contig master/slave coinvolti nel merge
		for( std::list<PairedContig>::const_iterator pctg = result.first.begin(); pctg != result.first.end(); pctg++ )
		{
			m_merged = (pctg->getMasterCtgMap()).size();
			s_merged = (pctg->getSlaveCtgMap()).size();

			m_merged_list.push_back( m_merged );
			s_merged_list.push_back( s_merged );

			m_merged_tot += m_merged;
			s_merged_tot += s_merged;

			merged_len += pctg->size();
			merged_num++;
		}

		_g_statsFile << "Master merged (mean) = " << double( double(m_merged_tot) / result.first.size() ) << "\n";
		_g_statsFile << " Slave merged (mean) = " << double( double(s_merged_tot) / result.first.size() ) << "\n";
		_g_statsFile << " Merged Total Length = " << merged_len << "\tAverage = " << merged_len / double(merged_num) << "\n";

		/*_g_statsFile << "Master merged list = ";
		for( std::list<uint32_t>::const_iterator m = m_merged_list.begin(); m != m_merged_list.end(); m++ )
			if( m == m_merged_list.begin() ) _g_statsFile << *m; else _g_statsFile << "\t" << *m;
		_g_statsFile << std::endl;

		_g_statsFile << "Slave merged list = ";
		for( std::list<uint32_t>::const_iterator m = s_merged_list.begin(); m != s_merged_list.end(); m++ )
			if( m == s_merged_list.begin() ) _g_statsFile << *m; else _g_statsFile << "\t" << *m;
		_g_statsFile << std::endl;*/

		_g_statsFile << std::endl;
	}

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
	std::cout << "[writing paired contigs on file: " << ( options.outputFilePrefix + ".fasta" ) << "]" << std::endl;

	result.first.sort( orderPctgsByName );

	std::ofstream outFasta((options.outputFilePrefix + ".fasta").c_str(),std::ios::out);
	std::list< PairedContig >::const_iterator i;
	for( i = result.first.begin(); i != result.first.end(); i++ ) outFasta << *i << std::endl;
	outFasta.close();

	// save paired contigs descriptors to file
	std::cout << "[writing paired contigs descriptors on file: " << ( options.outputFilePrefix + ".pctgs" ) << "]" << std::endl;
	std::fstream pctgDescFile( (options.outputFilePrefix + ".pctgs").c_str(), std::fstream::out );
	writePctgDescriptors( pctgDescFile, result.first, mcRef, scRef );
	pctgDescFile.close();

	// save IDs of (slave) contigs NOT merged
	std::cout << "[writing unused slave contigs on file: " << ( options.outputFilePrefix + ".unused" ) << "]" << std::endl;
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
	_g_statsFile.close();

	std::cout << "[total execution time = " << formatTime( time(NULL)-tStart ) << "]" << std::endl;

	std::ofstream ofs_fpi( (options.outputFilePrefix + ".fpi").c_str() );
	int max_frame_num = 0;
	for( int i=1; i < ext_fpi.size(); i++ ){ if( (ext_fpi[i] > 0 || ext_fpi_2[i] > 0) && i > max_frame_num) max_frame_num = i; }
	ofs_fpi << "          \t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << i << "\t"; } ofs_fpi << std::endl;
	ofs_fpi << "successful\t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << ext_fpi[i] << "\t"; } ofs_fpi << std::endl;
	ofs_fpi << "failed    \t";
	for( int i=1; i <= max_frame_num; i++ ){ ofs_fpi << ext_fpi_2[i] << "\t"; } ofs_fpi << std::endl;
	ofs_fpi.close();

	//g_badAlignStream.close();
}

} // namespace modules
