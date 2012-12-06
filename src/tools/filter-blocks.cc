/* 
 * File:   filter_blocks.cc
 * Author: Riccardo Vicedomini
 *
 * Created on 1 giugno 2011, 1.56
 */

#include<iostream>

#include <map>
#include <list>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <unistd.h>

#include <boost/graph/graphviz.hpp>

#include "types.hpp"
#include "graphs/PairingEvidencesGraph.hpp"
#include "pool/ContigMemPool.hpp"
#include "PartitionFunctions.hpp"

/*
 * 
 */
int main(int argc, char** argv) 
{
    if( argc != 5 )
    {
        std::cerr << "Usage: " << argv[0] << " <blocks file> <min pairing evidence> <master.fasta> <slave.fasta>" << std::endl;
        return 1;
    }
    
    
    std::vector< Block > blocks, outBlocks;
    std::string blockFile(argv[1]);
    int minEvidence = atoi(argv[2]);
    
    blocks = Block::readBlocks( blockFile );
    std::cerr << "Blocchi caricati: " << blocks.size() << std::endl;
    
    PairingEvidencesGraph peg(blocks);
    
    {
        std::ofstream f("PEG_test.dot");
        boost::write_graphviz(f,peg);
        f.close();
    }
    
    outBlocks = filterBlocksByPairingEvidences( blocks, minEvidence );
    std::cerr << "Blocchi filtrati: " << outBlocks.size() << std::endl;
    
    std::list< std::vector<Block> > pcblocks = partitionBlocks( outBlocks );
    std::cerr << "Paired Contigs: " << pcblocks.size() << std::endl;
    
    ContigMemPool masterPool, slavePool;
    
    masterPool = ContigMemPool::loadPool(argv[3]);
    slavePool = ContigMemPool::loadPool(argv[4]);
    
    
    
    return 0;
}

