/* 
 * File:   memory-usage.cc
 * Author: riccardo
 *
 * Created on 8 ottobre 2011, 20.53
 */

#include <iostream>
#include <sstream>
#include <fstream>

#include <string>
#include <map>

#include <time.h>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "assembly/io_contig.hpp"
#include "pool/ContigMemPool.hpp"
#include "pool/HashContigMemPool.hpp"

using namespace BamTools;

int main(int argc, char** argv) 
{    
    if( argc != 5 )
    {
        std::cout << "Usage: " << argv[0] << " <fasta1> <bam1> <fasta2> <bam2>" << std::endl;
        return 1;
    }
    
    time_t t1,t2;
    t1 = time(NULL);
    
    std::string fasta1 = argv[1];
    std::string bam1 = argv[2];
    std::string fasta2 = argv[3];
    std::string bam2 = argv[4];
    
    BamReader inBamM, inBamS;
    inBamM.Open( bam1 );
    inBamS.Open( bam2 );
    
    // load reference data (index => contig name)
    BamTools::RefVector mcRef = inBamM.GetReferenceData();
    BamTools::RefVector scRef = inBamS.GetReferenceData();

    //sparse_hash_map< std::string, int32_t > masterContigs( mcRef.size() ), slaveContigs( scRef.size() );
    //masterContigs.set_empty_key(std::string(""));
    //slaveContigs.set_empty_key(std::string(""));
    std::map< std::string, int32_t > masterContigs, slaveContigs;

    BamTools::RefVector::const_iterator ref_iter;

    for( ref_iter = mcRef.begin(); ref_iter != mcRef.end(); ref_iter++ ) 
        masterContigs[ ref_iter->RefName ] = ref_iter->RefLength;

    for( ref_iter = scRef.begin(); ref_iter != scRef.end(); ref_iter++ ) 
        slaveContigs[ ref_iter->RefName ] = ref_iter->RefLength;   
    
    std::cout << "Reference Data loaded." << std::endl;
    
    //USO DELLA STRUTTURA SPARSE-HASH
    HashContigMemPool pool1, pool2;
    //load master and slave contigs in memory
    pool1.loadPool(fasta1,masterContigs);
    std::cout << "Fasta1 loaded." << std::endl;
    pool2.loadPool(fasta2,slaveContigs);
    std::cout << "Fasta2 loaded." << std::endl;
    
    t2 = time(NULL);
    std::cout << "[Execution time: " << (t2-t1) << "]" << std::endl << std::flush;
    
    std::cout << "Press Enter to terminate..." << std::endl << std::flush;
    char c = getchar();

    return 0;
}

