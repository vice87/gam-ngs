/*
 * File:   hash_test.cc
 * Author: vice
 *
 * Created on March 26, 2012, 6:33 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include <google/dense_hash_set>

#include "UtilityFunctions.hpp"

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "bam/MultiBamReader.hpp"

using namespace BamTools;
using google::dense_hash_set;

int main(int argc, char** argv)
{
	if( argc != 4 )
	{
		std::cerr << "Usage: " << argv[0] << " <masterBAMs_list.txt> <slaveBAMs_list> <output.file>" << std::endl;
		return 1;
	}

	time_t start = time(NULL);

	int32_t slave_stats_id = -1, last_saved_ctg = -1;
	uint64_t slave_U = 0; // uniquely mapped reads on slaveCtgId
	uint64_t slave_M = 0; // multiple mapped reads on slaveCtgId
	uint64_t slave_U_master_N = 0; // uniquely mapped reads on slaveCtgId that are unmapped in the master
	uint64_t slave_U_master_U = 0; // uniquely mapped reads on slaveCtgId that are uniquely mapped in the master
	uint64_t slave_U_master_M = 0; // uniquely mapped reads on slaveCtgId that are multiple mapped in the master
	uint64_t slave_M_master_N = 0; // multiple mapped reads on slaveCtgId that are unmapped in the master
	uint64_t slave_M_master_U = 0;
	uint64_t slave_M_master_M = 0;
	uint64_t reads_len = 0;

	uint64_t g_slave_U = 0, g_slave_M = 0;
	uint64_t g_slave_U_master_N = 0, g_slave_U_master_U = 0, g_slave_U_master_M = 0;
	uint64_t g_slave_M_master_N = 0, g_slave_M_master_U = 0, g_slave_M_master_M = 0;

	uint64_t g_master_U = 0, g_master_M = 0;

	std::vector< int32_t > minInsert, maxInsert;

	uint64_t unmapped_contigs = 0, unmapped_contigs_len = 0;

	std::string output_file = argv[3];
	std::ofstream ofs( output_file.c_str() );

	MultiBamReader masterBamReader, slaveBamReader;
	std::vector< std::string > bamFiles; // vector of BAM filenames
	std::string filename;

	filename = argv[1];
	loadBamFileNames( filename, bamFiles, minInsert, maxInsert );
	masterBamReader.Open( bamFiles );

	filename = argv[2];
	bamFiles.clear();
	loadBamFileNames( filename, bamFiles, minInsert, maxInsert );
	slaveBamReader.Open( bamFiles );

	dense_hash_set< std::string >::const_iterator it;
	dense_hash_set< std::string > readSet_1, readSet_2;
	dense_hash_set< std::string > dupReadSet_1, dupReadSet_2;

	readSet_1.set_empty_key(std::string(""));
	readSet_2.set_empty_key(std::string(""));
	dupReadSet_1.set_empty_key(std::string(""));
	dupReadSet_2.set_empty_key(std::string(""));

	std::cout << "[loading master reads]" << std::endl;

	int32_t nh, xt; // molteplicità delle read (nh->standard, xt->bwa)

	BamAlignment align;

	while( masterBamReader.GetNextAlignment(align) )
	{
		// discard unmapped reads and reads that have a bad quality
		if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

		// load string fileds from bam
		// if( !align.BuildCharData() ) continue;

		// se la molteplicità non è stata definita, assumo che sia pari ad 1
		if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
		if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
		bool uniqMapRead = (nh == 1 && xt == 'U');

		if( uniqMapRead ) // read mappata in modo univoco
		{
			g_master_U++;
			// insert reads one of the reads hash-tables depending on whether it is the first or second pair
			if( align.IsFirstMate() || !align.IsPaired() ) readSet_1.insert(align.Name); else readSet_2.insert(align.Name);
		}
		else
		{
			g_master_M++;
			// insert reads one of the reads hash-tables depending on whether it is the first or second pair
			if( align.IsFirstMate() || !align.IsPaired() ) dupReadSet_1.insert(align.Name); else dupReadSet_2.insert(align.Name);
		}
	}

	masterBamReader.Close();
	std::cout << "[reads loaded in " << formatTime( time(NULL)-start ) << "]" << std::endl;
	std::cout << "[computing statistics]" << std::endl;

	const RefVector& refVect = slaveBamReader.GetReferenceData();

	while( slaveBamReader.GetNextAlignment(align) )
	{
		if( align.IsMapped() && align.RefID < slave_stats_id && slave_stats_id != -1 )
		{
			std::cerr << "Alignments aren't being retrived by contig order. This should not happend if BAMs are sorted." << std::endl;
			return 1;
		}

		if( align.IsMapped() && align.RefID > slave_stats_id ) // moving to a new contig
		{
			// save statistics (only if id is valid)
			if( slave_stats_id != -1 )
			{
				ofs << refVect[slave_stats_id].RefName << "\t" << slave_U << "\t" << slave_M << "\t"
				<< slave_U_master_N << "\t" << slave_U_master_U << "\t" << slave_U_master_M << "\t"
				<< slave_M_master_N << "\t" << slave_M_master_U << "\t" << slave_M_master_M << "\t"
				<< refVect[slave_stats_id].RefLength << std::endl;

				if( double(slave_U_master_N) >= 0.9 * double(slave_U) )
				{
					unmapped_contigs++;
					unmapped_contigs_len += refVect[slave_stats_id].RefLength;

					std::cout << refVect[slave_stats_id].RefName << "\t" << slave_U << "\t" << slave_M << "\t"
					<< slave_U_master_N << "\t" << slave_U_master_U << "\t" << slave_U_master_M << "\t"
					<< slave_M_master_N << "\t" << slave_M_master_U << "\t" << slave_M_master_M << "\t"
					<< refVect[slave_stats_id].RefLength << std::endl;
				}

				last_saved_ctg = slave_stats_id;
			}

			// update counters for new contig
			slave_stats_id = align.RefID;
			slave_U = slave_M = 0;
			slave_U_master_N = slave_U_master_U = slave_U_master_M = 0;
			slave_M_master_N = slave_M_master_U = slave_M_master_M = 0;
			reads_len = 0;
		}

		if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

		// se la molteplicità non è stata definita, assumo che sia pari ad 1
		if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard SAM format field
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // bwa field
		bool uniqMapRead = (nh == 1 && xt == 'U');

		if( uniqMapRead )
		{
			slave_U++; g_slave_U++;
			reads_len += (align.GetEndPosition() - align.Position);

			if( !align.IsPaired() || align.IsFirstMate() ) // cerca la read corrente fra quelle del master
			{
				it = readSet_1.find(align.Name);
				if( it != readSet_1.end() )
				{
					slave_U_master_U++; g_slave_U_master_U++;
				}
				else
				{
					it = dupReadSet_1.find(align.Name);
					if( it != dupReadSet_1.end() ){ slave_U_master_M++; g_slave_U_master_M++; } else { slave_U_master_N++; g_slave_U_master_N++; }
				}
			}
			else
			{
				it = readSet_2.find(align.Name);
				if( it != readSet_2.end() )
				{
					slave_U_master_U++; g_slave_U_master_U++;
				}
				else
				{
					it = dupReadSet_2.find(align.Name);
					if( it != dupReadSet_2.end() ){ slave_U_master_M++; g_slave_U_master_M++; } else { slave_U_master_N++; g_slave_U_master_N++; }
				}
			}
		}
		else // read not uniquely mapped
		{
			slave_M++; g_slave_M++;

			if( !align.IsPaired() || align.IsFirstMate() ) // cerca la read corrente fra quelle del master
			{
				it = readSet_1.find(align.Name);
				if( it != readSet_1.end() )
				{
					slave_M_master_U++; g_slave_M_master_U++;
				}
				else
				{
					it = dupReadSet_1.find(align.Name);
					if( it != dupReadSet_1.end() ){ slave_M_master_M++; g_slave_M_master_M++; } else { slave_M_master_N++; g_slave_M_master_N++; }
				}
			}
			else
			{
				it = readSet_2.find(align.Name);
				if( it != readSet_2.end() )
				{
					slave_M_master_U++; g_slave_M_master_U++;
				}
				else
				{
					it = dupReadSet_2.find(align.Name);
					if( it != dupReadSet_2.end() ){ slave_M_master_M++; g_slave_M_master_M++; } else { slave_M_master_N++; g_slave_M_master_N++; }
				}
			}
		}

		/////
	}

	// save statistics for last contig
	if( slave_stats_id != -1 && slave_stats_id != last_saved_ctg )
	{
		ofs << refVect[slave_stats_id].RefName << "\t" << slave_U << "\t" << slave_M << "\t"
			<< slave_U_master_N << "\t" << slave_U_master_U << "\t" << slave_U_master_M << "\t"
			<< slave_M_master_N << "\t" << slave_M_master_U << "\t" << slave_M_master_M << "\t"
			<< refVect[slave_stats_id].RefLength << std::endl;

		if( double(slave_U_master_N) >= 0.9 * double(slave_U) )
		{
			unmapped_contigs++;
			unmapped_contigs_len += refVect[slave_stats_id].RefLength;

			std::cout << refVect[slave_stats_id].RefName << "\t" << slave_U << "\t" << slave_M << "\t"
			<< slave_U_master_N << "\t" << slave_U_master_U << "\t" << slave_U_master_M << "\t"
			<< slave_M_master_N << "\t" << slave_M_master_U << "\t" << slave_M_master_M << "\t"
			<< refVect[slave_stats_id].RefLength << std::endl;
		}
	}

	ofs.close();

	slaveBamReader.Close();

	time_t end = time(NULL);

	std::cout << "[execution time = " << formatTime( end-start ) << "]" << std::endl;

	std::cout << "\n[master reads'mapping statistics:]" << std::endl;
	std::cout << "unique = " << g_master_U << "\tmultiple = " << g_master_M << std::endl;

	std::cout << "\n[slave statistics:]" << std::endl;
	std::cout << "unique = " << g_slave_U << "\tmultiple = " << g_slave_M << std::endl;
	std::cout << "unique => unique = " << g_slave_U_master_U << "\tmultiple = " << g_slave_U_master_M << "\tunmapped = " << g_slave_U_master_N << std::endl;
	std::cout << "multiple => unique = " << g_slave_M_master_U << "\tmultiple = " << g_slave_M_master_M << "\tunmapped = " << g_slave_M_master_N << std::endl;

	std::cout << "\nunmapped contigs = " << unmapped_contigs << "\tlength = " << unmapped_contigs_len << std::endl;

    return 0;
}

