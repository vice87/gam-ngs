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
    uint32_t reads = 0, mapped = 0, uniMapped = 0, peCorrectlyMapped = 0;
    int8_t type;
        
    while( inBam.GetNextAlignment(align) )
    {
        reads++;
        
        if( align.IsMapped() ) mapped++;
        
        // uniquely mapped reads
        if( align.IsMapped() && align.GetTag( std::string("XT"), type ) )
        {
            if( type == 'U' ) uniMapped++;
        }
        
        if( align.IsMapped() )
        {
            if( align.IsPaired() && align.IsProperPair() && align.InsertSize >= minTemplateLen && align.InsertSize <= maxTemplateLen ) peCorrectlyMapped++;
        }
    }
    
    inBam.Close();
    
    std::cout << "Reads: " << reads << "\n"
              << "Reads Mapped: " << mapped << "\n"
              << "Reads Uniquely Mapped: " << uniMapped << "\n"
              << "PE Correctly Mapped: " << peCorrectlyMapped << "\n"
              << std::endl;
    
    return 0;
}

