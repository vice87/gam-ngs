/* 
 * File:   unused-contigs.cc
 * Author: vice
 *
 * Created on 10 novembre 2011, 16.22
 */

#include <cstdlib>
#include "api/BamAux.h"
#include "api/BamReader.h"
#include "types.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) 
{
    if( argc != 3 )
    {
        std::cout << "Usage: " << argv[0] << " [pctgs.file] [slave.bam.file] > [output.file]" << std::endl;
        return 1;
    }
    
    std::string pctgsFile = argv[1];
    std::string bamFile = argv[2];
    
    BamTools::BamReader inBam;
    inBam.Open( bamFile );
    
    BamTools::RefVector bamSequences = inBam.GetReferenceData();
    std::vector<bool> usedCtgs( bamSequences.size(), false );
    
    std::ifstream input( pctgsFile.c_str() );
    
    while( !input.eof() )
    { 
        char c = input.peek();
        
        while( c == '\n' || c == '\r' || c == '\t' ){ input.ignore(1); c = input.peek(); }
        
        if( c == '#' )
        {
            while( c != '\n' ){ input.ignore(1); c = input.peek(); }
            input.ignore(1);
            continue;
        }
        else
        {
            std::string trash;
            std::string type;
            UIntType ctgID;
            input >> trash >> trash >> type >> ctgID >> trash >> trash;
            
            if( type == "Slave" ) usedCtgs[ctgID] = true;
        }
    }
    
    input.close();
    
    for( UIntType i = 0; i < usedCtgs.size(); i++ )
    {
        if( !usedCtgs[i] ) std::cout << bamSequences[i].RefName << std::endl;
    }
    
    return 0;
}

