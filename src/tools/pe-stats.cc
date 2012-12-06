/*
 * File:   pe-stats.cc
 * Author: riccardo
 *
 * Created on 1 agosto 2011, 22.39
 */

#include <cstdlib>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>

#include "api/BamReader.h"

using namespace BamTools;

std::string getPathBaseName( std::string path )
{
    size_t found = path.rfind ( '/' );
    return (found != std::string::npos) ? path.substr(found) : std::string( path.c_str() );
}

int main(int argc, char** argv)
{
    if( argc != 3 )
    {
        std::cout << "Usage " << getPathBaseName( std::string(argv[0]) ) << "<insert length> <alignment.bam>" << std::endl;
        return 1;
    }

    std::string bamFile = argv[2];
    int32_t insertLength = atoi(argv[1]);

    int32_t minTemplateLen = insertLength - (insertLength/2);
    int32_t maxTemplateLen = insertLength + (insertLength/2);

    std::cout << "Insert Size = " << minTemplateLen << " - " << maxTemplateLen << std::endl;

    BamReader inBam;
    inBam.Open( bamFile );

    BamAlignment align;
    uint64_t reads = 0, mapped = 0, mappedQC = 0, mappedQC_uniq = 0, mappedQC_uniq_paired = 0, mappedQc_uniq_paired_correctlyMapped = 0;
    int8_t type;
	int32_t nh,xt;

    while( inBam.GetNextAlignment(align) )
    {
        reads++;

		if( !align.IsMapped() ) continue;
		mapped++;

		if( align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;
		mappedQC++;

		if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
		if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
		if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1
		mappedQC_uniq++;

		if( !align.IsPaired() ) continue;
		mappedQC_uniq_paired++;

		if( !align.IsProperPair() || align.InsertSize < minTemplateLen || align.InsertSize > maxTemplateLen ) continue;
		mappedQc_uniq_paired_correctlyMapped++;
    }

    inBam.Close();

    std::cout << "Reads: " << reads << "\n"
              << "Reads Mapped: " << mapped << "\n"
              << "Reads Mapped (good quality): " << mappedQC << "\n"
			  << "Reads Mapped (good quality), unique: " << mappedQC_uniq << "\n"
			  << "Reads Mapped (good quality), unique, paired: " << mappedQC_uniq_paired << "\n"
			  << "Reads Mapped (good quality), unique, paired, correctly mapped: " << mappedQc_uniq_paired_correctlyMapped*2 << "\n"
              << std::endl;

    return 0;
}

