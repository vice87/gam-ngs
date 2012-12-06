/* 
 * File:   hash_test.cc
 * Author: vice
 *
 * Created on March 26, 2012, 6:33 PM
 */

#include <cstdlib>
#include <iostream>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using namespace BamTools;


bool included_in( BamAlignment &align, int32_t start, int32_t end )
{
    int32_t read_start = align.Position;
    int32_t read_end = align.GetEndPosition() - 1;
    
    return( start <= read_start && read_end <= end );
}


int main(int argc, char** argv) 
{
    if( argc != 7 )
    {
        std::cout << "Usage: pe-test <input.bam> <ctg_id> <s1> <e1> <s2> <e2>" << std::endl;
        return 1;
    }
    
    std::string inputBam = argv[1];
    int32_t id = atoi( argv[2] );
    int32_t s1 = atoi( argv[3] );
    int32_t e1 = atoi( argv[4] );
    int32_t s2 = atoi( argv[5] );
    int32_t e2 = atoi( argv[6] );
    
    if( s1 > e1 || s2 > e2 )
    {
        std::cout << "Wrong interval definition. It has to be s1 <= e1 and s2 <= e2" << std::endl;
        return 1;
    }
    
    BamReader bamReader;
    BamAlignment align;
    
    if( !bamReader.Open( inputBam ) ){ std::cout << "cannot open bam file" << std::endl; return 1; }
    if( !bamReader.OpenIndex( inputBam + ".bai" ) ){ std::cout << "cannot open bam index file" << std::endl; return 1; }

    int32_t nh, xt;
    
    uint64_t collisions = 0;
    uint64_t num_reads = 0;
    uint64_t mates = 0;
    
    bamReader.SetRegion( id, s1, id, e2 );
    
    while( bamReader.GetNextAlignmentCore(align) )
    {
        if( !align.IsPaired() ) continue;
        if( !align.IsMapped() || !align.IsMateMapped() ) continue;
        
        if( align.RefID != align.MateRefID ) continue; // pairs must align in the same contig
        
        align.BuildCharData(); // fill string fields
        
        // se la molteplicità non è stata definita, assumo che sia pari ad 1
        if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
        
        if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1
        
        if( included_in(align,s1,e1) && align.MatePosition >= s2 && align.MatePosition <= e2 ) mates++; 
        if( included_in(align,s2,e2) && align.MatePosition >= s1 && align.MatePosition <= e1 ) mates++;
    }
    
    bamReader.Close();
    
    std::cout << "Number of pairs that bridge the blocks: " << mates << std::endl;
    
    return 0;
}

