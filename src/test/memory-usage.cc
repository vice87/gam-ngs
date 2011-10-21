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
    
    std::string fasta1 = argv[1];
    std::string bam1 = argv[2];
    std::string fasta2 = argv[3];
    std::string bam2 = argv[4];
    
    //USO DELLA STRUTTURA SPARSE-HASH
    HashContigMemPool pool1, pool2;
    //load master and slave contigs in memory
    pool1.loadPool(fasta1);
    std::cout << "Fasta1 loaded." << std::endl;
    pool2.loadPool(fasta2);
    std::cout << "Fasta2 loaded." << std::endl;
    
    //USO DI UN SEMPLICE ARRAY
//    BamReader inBamM, inBamS;
//    inBamM.Open( bam1 );
//    inBamS.Open( bam2 );
//    
//    BamTools::RefVector mcRef = inBamM.GetReferenceData();
//    BamTools::RefVector scRef = inBamS.GetReferenceData();
//    
//    Contig ctgs1[ mcRef.size() ];
//    Contig ctgs2[ scRef.size() ];
//    
//    std::ifstream ifs( fasta1.c_str(), std::ifstream::in );
//    long count = 0;
//    
//    while( !ifs.eof() )
//    {
//        Contig ctg;
//        
//        ifs >> ctg;
//        ctgs1[count] = ctg;
//        count++;
//    }
//    
//    ifs.open( fasta2.c_str(), std::ifstream::in );
//    
//    count = 0;
//    
//    while( !ifs.eof() )
//    {
//        Contig ctg;
//        
//        ifs >> ctg;
//        ctgs2[count] = ctg;
//        count++;
//    }
//    
//    ifs.close();
    
    std::cout << "Press Enter to terminate..." << std::endl << std::flush;
    char c = getchar();

    return 0;
}

