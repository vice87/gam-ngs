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
	time_t start = time(NULL);

	std::string filename = argv[1];
	std::vector< std::string > bamFiles; // vector of BAM filenames
	std::vector< int32_t > minInsert, maxInsert;

	MultiBamReader bamReader;
	loadBamFileNames( filename, bamFiles, minInsert, maxInsert );
	bamReader.Open( bamFiles ); // open master BAM files

	const RefVector& ref = bamReader.GetReferenceData();

	std::vector< uint64_t > reads_len( ref.size(), 0 );
	std::vector< uint64_t > coverage( 51, 0 );

	std::vector< uint64_t > reads_len_2( ref.size(), 0 );
	std::vector< uint64_t > coverage_2( 51, 0 );

	int32_t nh, xt; // molteplicità delle read (nh->standard, xt->bwa)

	BamAlignment align;
	while( bamReader.GetNextAlignment(align) )
	{
		if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

		int32_t length = align.GetEndPosition() - align.Position;

		reads_len_2[ align.RefID ] += length;

		// se la molteplicità non è stata definita, assumo che sia pari ad 1
		if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard SAM format field
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // bwa field

        // scarto read con molteplicità maggiore di 1
		if( nh != 1 || xt != 'U' ) continue;

		reads_len[ align.RefID ] += length;
	}

	// compute coverage for each reference sequence
	for( size_t i=0; i < ref.size(); i++ )
	{
		double c = reads_len[i] / ref[i].RefLength;
		uint32_t idx = (c < 0) ? 0 : uint32_t(c/2.0);
		if( idx < 50 ) coverage[idx]++; else coverage[50]++;

		double c2 = reads_len_2[i] / ref[i].RefLength;
		uint32_t idx2 = (c2 < 0) ? 0 : uint32_t(c2/2.0);
		if( idx2 < 50 ) coverage_2[idx2]++; else coverage_2[50]++;
	}

	std::string output_file = filename + ".coverage.stats";
	std::ofstream ofs( output_file.c_str() );
	for( size_t i=0; i < coverage.size(); i++ ) ofs << i*2 << "\t" << coverage[i] << "\n";
	ofs.close();

	output_file = filename + ".coverage_nu.stats";
	ofs.open( output_file.c_str() );
	for( size_t i=0; i < coverage_2.size(); i++ ) ofs << i*2 << "\t" << coverage_2[i] << "\n";
	ofs.close();

	bamReader.Close();

	time_t end = time(NULL);

	std::cout << "time = " << formatTime( end-start ) << std::endl;

    return 0;
}

