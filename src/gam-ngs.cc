/* 
 * File:   samgam.cc
 * Author: RiccardoVicedomini
 *
 * Created on 21 maggio 2011, 14.53
 */

#include<iostream>

#include <map>
#include <list>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

#include <boost/graph/graphviz.hpp>
#include <google/sparse_hash_map>

#include "types.hpp"
#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "assembly/Read.hpp"
#include "assembly/Block.hpp"
#include "graphs/PairingEvidencesGraph.hpp"
#include "pctg/PairedContig.hpp"
#include "pool/HashContigMemPool.hpp"
#include "OrderingFunctions.hpp"
#include "PartitionFunctions.hpp"
#include "ThreadedBuildPctg.code.hpp"

using namespace BamTools;
using google::sparse_hash_map;

/*
 * 
 */
int main(int argc, char** argv) 
{
    int minBlockSize;           // dimensione minima per i blocchi
    std::string bamFileM;	// file BAM (Master)
    //std::string bamIndexFileM;	// file indice del file BAM (Master)
    std::string bamMasterSorted;
    std::string bamFileS;	// file BAM (Slave)
    //std::string bamIndexFileS;	// file indice del file BAM (Slave)
    std::string bamSlaveSorted;
    std::string masterFasta;
    std::string slaveFasta;
    std::string action;
    int minEvidence;
    size_t threadsNum;
    
    switch(argc)
    {
        case 11:
            threadsNum = atoi(argv[10]);
            if( threadsNum <= 0 ) threadsNum = 1;
        case 10: 
            action = argv[1];
            bamFileM = argv[2];
            //bamIndexFileM = argv[2];
            bamMasterSorted = argv[3];
            bamFileS = argv[4];
            //bamIndexFileS = argv[5];
            bamSlaveSorted = argv[5];
            minBlockSize = atoi(argv[6]);
            minEvidence = atoi(argv[7]);
            masterFasta = argv[8];
            slaveFasta = argv[9];
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " [blocks|merge]"
                      << "<BAM Master CoordinateSorted> <BAM Master NameSorted> "
                      << "<BAM Slave CoordinateSorted> <BAM Slave NameSorted> "
                      << "<Min Block Size> <Min Pairing Evidence> "
                      << "<Master.fasta> <Slave.fasta> <Number of Threads>"
                      << std::endl;
            return 1;
            break;
    }    
    
    if( action.compare("blocks") != 0 && action.compare("merge") != 0 )
    {
        std::cerr << "Usage: " << argv[0] << " [blocks|merge]"
                  << "<BAM Master CoordinateSorted> <BAM Master NameSorted> "
                  << "<BAM Slave CoordinateSorted> <BAM Slave NameSorted> "
                  << "<Min Block Size> <Min Pairing Evidence> "
                  << "<Master.fasta> <Slave.fasta> <Number of Threads>"
                  << std::endl;
        
        return 1;
    }
    
    time_t t1,t2;
    t1 = time(NULL);
    
    BamReader inBamM;
    inBamM.Open( bamFileM );
    
    struct stat st;
    if( stat((bamFileM + ".bai").c_str(),&st) == 0 ) inBamM.OpenIndex( bamFileM + ".bai" );
    
    BamReader inBamS;
    inBamS.Open( bamFileS );
    if( stat((bamFileS + ".bai").c_str(),&st) == 0 ) inBamM.OpenIndex( bamFileS + ".bai" );
    
    std::vector<Block> blocks;
    
    if( action.compare("blocks") == 0 ) // create blocks file
    {
        BamReader inBamMasterSorted;
        inBamMasterSorted.Open( bamMasterSorted );

        BamReader inBamSlaveSorted;
        inBamSlaveSorted.Open( bamSlaveSorted );
        
        // load only useful reads of the slave assembly
        sparse_hash_map< std::string, Read > readMap;
        Read::getReadMap(inBamMasterSorted,inBamSlaveSorted,readMap);
        
        inBamMasterSorted.Close();
        inBamSlaveSorted.Close();
        
        // find and build blocks.
        blocks = Block::findBlocks( inBamM, inBamS, minBlockSize, readMap );
        std::cerr << "Number of blocks found: " << blocks.size() << std::endl;
        
        // write found blocks to file
        Block::writeBlocks( std::string("blocks.txt"), blocks );
    }
    
    if( action.compare("merge") == 0 ) // load blocks from file
    {
        blocks = Block::readBlocks( std::string("blocks.txt") );
        std::cerr << "Blocchi caricati: " << blocks.size() << std::endl;
    }
    
    //PairingEvidencesGraph peg(blocks);
    //
    //{
    //    std::ofstream f("PEG_test.dot");
    //    boost::write_graphviz(f,peg);
    //    f.close();
    //}
    
    std::vector< Block > outBlocks;
    
    // keep only blocks between contigs that share at least minEvidence blocks.
    outBlocks = filterBlocksByPairingEvidences( blocks, minEvidence );
    std::cerr << "Filtered Blocks: " << outBlocks.size() << std::endl << std::flush;
    
    // create the graph of assemblies, remove cycles and keep remaining blocks.
    std::list< std::vector<Block> > pcblocks = partitionBlocks( outBlocks );
    std::cerr << "Paired Contigs Blocks: " << pcblocks.size() << std::endl << std::flush;
    
    std::cerr << "Loading contigs in memory... " << std::flush;
    
    HashContigMemPool masterPool, slavePool, pctgPool;
    
    // load master and slave contigs in memory
    masterPool.loadPool(masterFasta);
    slavePool.loadPool(slaveFasta);
    
    std::cerr << "done." << std::endl << std::flush;
    
    BamTools::RefVector mcRef = inBamM.GetReferenceData();
    BamTools::RefVector scRef = inBamS.GetReferenceData();
    
    std::cerr << "ReferenceData loaded: " << mcRef.size() << " master contigs, " << scRef.size() << " slave contigs" << std::endl << std::flush;
    
    // build paired contigs (threaded)
    ThreadedBuildPctg tbp( pcblocks, &pctgPool, &masterPool, &slavePool, &mcRef, &scRef );
    std::pair< std::list<PairedContig>, std::vector<bool> > result = tbp.run(threadsNum);
    
    IdType pctgNum((result.first).size());
    std::cerr << "Paired Contigs: " << pctgNum << std::endl << std::flush;
    
    // slave contig pool is no more needed
    slavePool.clear();
    
    std::list<IdType> ctgIds;
    
    for(IdType i = 0; i < result.second.size(); i++) 
        if( !result.second.at(i) ) ctgIds.push_back(i);
    
    generateSingleCtgPctgs( result.first, ctgIds, &pctgPool, &masterPool, &mcRef, pctgNum);
    
    // save paired contig pool to file
    pctgPool.savePool(std::string("pctgs.fasta"));
    
    // save paired contigs descriptors to file
    std::cerr << "Writing PairedContig descriptors..." << std::flush;
    
    result.first.sort( orderPctgsByName );
    
    std::fstream pctgDescFile("pctgs-descriptors.txt",std::fstream::out);
    writePctgDescriptors( pctgDescFile, result.first );
    pctgDescFile.close();
    
    std::cerr << "done" << std::endl;
    
    t2 = time(NULL);
    
    std::cerr << "Execution time: " << t2-t1 << "s" << std::endl;
    
    inBamM.Close();
    inBamS.Close();
    
    return 0;
}
