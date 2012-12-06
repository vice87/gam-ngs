/* 
 * File:   hash_test.cc
 * Author: vice
 *
 * Created on March 26, 2012, 6:33 PM
 */

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "murmur3.hpp"
#include <google/sparse_hash_map>

#include "assembly/Read.hpp"
#include "UtilityFunctions.hpp"

using google::sparse_hash_map;

int main(int argc, char** argv) 
{
    std::string read_name;
    std::size_t len;
    std::size_t hash[2];
    char line[1024];
    char *ptr;
    Read a_read;
    std::pair< std::size_t, Read > pair;
    
    uint64_t collisions = 0;
    uint64_t num_reads = 0;
    
    std::cout << "Size of double = " << sizeof(double) << std::endl;
    
    google::sparse_hash_map< std::size_t, Read > readMap_1, readMap_2;
    
    for( int i=1; i < argc; i++ )
    {
        std::ifstream input( argv[i] );
        std::cout << "Processing " << getPathBaseName( argv[i] ) << std::endl;
        
        while( input.good() )
        {
            input.getline( line, 1024 );
            
            read_name = line;
            
            if( read_name.length() == 0 ) continue;
            if( read_name[0] != '@' ) continue;
            
            std::size_t p = read_name.find(' ');
            read_name = read_name.substr(0,p);
	    len = read_name.length();
            
            if( read_name.substr(len-2) == "/1" )
            {
                read_name = read_name.substr(0,len-2);
                MurmurHash3_x64_128( (void *) read_name.c_str(), (int) len-2, (uint32_t) 9001, (void *) hash  );
                pair = std::make_pair( hash[0], a_read );
                
                if( readMap_1.find( hash[0] ) == readMap_1.end() ) readMap_1.insert(pair); else collisions++;
            }
            else if( read_name.substr(len-2) == "/2" )
            {
                read_name = read_name.substr(0,len-2);
                MurmurHash3_x64_128( (void *) read_name.c_str(), (int) len-2, (uint32_t) 9001, (void *) hash  );
                pair = std::make_pair( hash[0], a_read );
                
                if( readMap_2.find( hash[0] ) == readMap_2.end() ) readMap_2.insert(pair); else collisions++;
            }
            
            // discard next 3 lines
            if( input.good() ) input.getline( line, 1024 );
            if( input.good() ) input.getline( line, 1024 );
            if( input.good() ) input.getline( line, 1024 );
        }
        
        std::cout << getPathBaseName( argv[i] ) << " loaded successfully!" << std::endl;
        
        double ratio = double(collisions) / double(num_reads);
        std::cout << "Total reads (uniquely mapped): " << num_reads << std::endl;
        std::cout << "Collisions: " << collisions << std::endl;
        std::cout << "Ratio: " << ratio << std::endl;
        
        std::cout << "Press Enter to continue..."  << std::endl;
        char c = getchar();
        
        input.close();
    }
      
    double ratio = double(collisions) / double(num_reads);
    
    std::cout << "Total reads (uniquely mapped): " << num_reads << std::endl;
    std::cout << "Collisions: " << collisions << std::endl;
    std::cout << "Ratio: " << ratio << std::endl;
    
    return 0;
}

