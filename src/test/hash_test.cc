/* 
 * File:   hash_test.cc
 * Author: vice
 *
 * Created on March 26, 2012, 6:33 PM
 */

#include <cstdlib>
#include "murmur3.hpp"
#include <google/sparse_hash_map>

#include "UtilityFunctions.hpp"
#include <time.h>

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "assembly/Read.hpp"

using namespace BamTools;
using google::sparse_hash_map;

int main(int argc, char** argv) 
{
    std::string inputBam = argv[1];
    
    BamReader bamReader;
    BamAlignment align;
    
    std::map< std::size_t, Read > readMap_1, readMap_2;
    //readMap_1.set_deleted_key("");
    //readMap_2.set_deleted_key("");
    
    size_t hash[2];
    
    time_t tStart = time(NULL);
    
    bamReader.Open( inputBam );
    int32_t nh, xt;
    
    uint64_t collisions = 0;
    uint64_t num_reads = 0;
    
    while( bamReader.GetNextAlignment(align) )
    {
        // se la molteplicità non è stata definita, assumo che sia pari ad 1
        if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
        
        if( nh == 1 && xt == 'U' && align.IsMapped() )
        {
            Read curRead( align.RefID, align.Position, align.GetEndPosition(), align.IsReverseStrand() );
            
            MurmurHash3_x64_128( (void *) align.Name.c_str(), (int) align.Name.length(), (uint32_t) 9001, (void *) hash  );
            std::pair< std::size_t, Read > pair = std::make_pair( hash[0], curRead );
            
            num_reads++;
            
            if( align.IsFirstMate() )
            {
                if( readMap_1.find( hash[0] ) == readMap_1.end() ) readMap_1.insert(pair); else collisions++;
            }
            else
            {
                if( readMap_2.find( hash[0] ) == readMap_2.end() ) readMap_2.insert(pair); else collisions++;
            }
        }
    }
    
    bamReader.Close();
    
    std::cout << "Execution Time: " << formatTime( time(NULL) - tStart ) << std::endl;
    
    double ratio = double(collisions) / double(num_reads);
    
    std::cout << "Total reads (uniquely mapped): " << num_reads << std::endl;
    std::cout << "Collisions: " << collisions << std::endl;
    std::cout << "Ratio: " << ratio << std::endl;
    
    return 0;
}

