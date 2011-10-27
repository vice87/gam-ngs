/* 
 * File:   samgam.cc
 * Author: RiccardoVicedomini
 *
 * Created on 21 maggio 2011, 14.53
 */

#include <iostream>

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
#include "UtilityFunctions.hpp"

using namespace BamTools;
using google::sparse_hash_map;

/*
 * 
 */
int main(int argc, char** argv) 
{
    std::string action;                 // command
    std::string outFilePrefix;          // output file prefix
    std::string inBlocksFile;           // input blocks file
    
    std::string bamFileM;               // file BAM (Master)
    std::string bamMasterSorted;        // file BAM (Master) sorted by read name
    std::string bamFileS;               // file BAM (Slave)
    std::string bamSlaveSorted;         // file BAM (Slave) sorted by read name
    
    int minBlockSize;                   // dimensione minima (# reads) per i blocchi
    int minEvidence;
    
    std::string masterFasta;            // FASTA file of the master assembly
    std::string slaveFasta;             // FASTA file of the slave assembly
    
    size_t threadsNum;                  // Numero di thread da usare per il weaving
    
    time_t t1,t2;
    t1 = time(NULL);
    
    if( argc < 2 )
    {
        std::cerr << std::endl << "Usage:   " << getPathBaseName(argv[0]) << " <command> [options]" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Command: " << "blocks" << "\t" << "generate blocks file" << std::endl
                  << "         " << "merge" << "\t" << "merge two assemblies given a blocks file" << std::endl
                  << std::endl;
        
        return 1;
    }
    
    action = argv[1]; // get command string
    
    if( action.compare("blocks") != 0 && action.compare("merge") != 0 )
    {
        std::cerr << std::endl << "Usage:   " << getPathBaseName(argv[0]) << " <command> [options]" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Command: " << "blocks" << "\t" << "generate blocks file" << std::endl
                  << "         " << "merge" << "\t" << "merge two assemblies given a blocks file" << std::endl
                  << std::endl;
        
        return 1;
    }
    
    /* Command to build the blocks */
    if( action.compare("blocks") == 0 )
    {
        if( argc != 8 )
        {
            std::cerr << "Usage: " << getPathBaseName(argv[0]) << " blocks "
                      << "<BAM Master CoordinateSorted> <BAM Master ReadNameSorted> " 
                      << "<BAM Slave CoordinateSorted> <BAM Slave ReadNameSorted> "
                      << "<Min Block Size> <Output FileName>"
                      << std::endl;
            return 1;
        }
        
        bamFileM = argv[2];
        bamMasterSorted = argv[3];
        bamFileS = argv[4];
        bamSlaveSorted = argv[5];
        minBlockSize = atoi(argv[6]);
        outFilePrefix = argv[7];
        
        if( minBlockSize < 1 ) minBlockSize = 1;
        
        BamReader inBamMasterSorted;
        inBamMasterSorted.Open( bamMasterSorted );

        BamReader inBamSlaveSorted;
        inBamSlaveSorted.Open( bamSlaveSorted );
        
        std::cerr << "[building reads map]" << std::endl;
        
        /* load only useful reads of the slave assembly */
        sparse_hash_map< std::string, Read > readMap;
        readMap.set_deleted_key("");
        Read::getReadMap(inBamMasterSorted,inBamSlaveSorted,readMap);
        
        inBamMasterSorted.Close();
        inBamSlaveSorted.Close();
        
        struct stat st;
        
        BamReader inBamM;
        inBamM.Open( bamFileM );
        if( stat((bamFileM + ".bai").c_str(),&st) == 0 ) inBamM.OpenIndex( bamFileM + ".bai" );

        BamReader inBamS;
        inBamS.Open( bamFileS );
        if( stat((bamFileS + ".bai").c_str(),&st) == 0 ) inBamM.OpenIndex( bamFileS + ".bai" );
        
        std::cerr << "[finding blocks]" << std::endl;
        std::vector<Block> blocks = Block::findBlocks( inBamM, inBamS, minBlockSize, readMap );
        
        std::cerr << "[writing blocks on file: " << getPathBaseName( outFilePrefix + ".blocks" ) << "]" << std::endl;
        Block::writeBlocks( outFilePrefix + ".blocks" , blocks ); // write blocks to file
        
        inBamM.Close();
        inBamS.Close();
        
        return 0;
    }
    
    /* Command to merge assemblies */
    if( action.compare("merge") == 0 )
    {
        if( argc != 10 )
        {
            std::cerr << "Usage: " << getPathBaseName(argv[0]) << " merge "
                      << "<Input.blocks> <Min Block Size> "
                      << "<BAM Master CoordinateSorted> <BAM Slave CoordinateSorted> "
                      << "<Master FASTA> <Slave FASTA> <Output FileName> "
                      << "<Threads>"
                      << std::endl;
            return 1;
        }
        
        inBlocksFile = argv[2];
        minBlockSize = atoi(argv[3]);
        bamFileM = argv[4];
        bamFileS = argv[5];
        masterFasta = argv[6];
        slaveFasta = argv[7];
        outFilePrefix = argv[8];
        threadsNum = atoi(argv[9]);
        
        if( threadsNum < 1 ) threadsNum = 1;
        
        BamReader inBamM, inBamS;
        inBamM.Open( bamFileM );
        inBamS.Open( bamFileS );
        
        std::cerr << "[loading blocks]" << std::endl;
        std::vector<Block> blocks = Block::readBlocks( inBlocksFile, minBlockSize );
        std::cerr << "[blocks loaded: " << blocks.size() << "]" << std::endl;
        
        // keep only blocks between contigs that share at least minEvidence blocks.
        // std::vector< Block > outBlocks = filterBlocksByPairingEvidences( blocks, minEvidence );
        
        // create the graph of assemblies, remove cycles and keep remaining blocks.
        std::cerr << "[partitioning blocks for merging]" << std::endl;
        std::list< std::vector<Block> > pcblocks = partitionBlocks( blocks );
        
        std::cerr << "[loading contigs in memory]" << std::endl;
        HashContigMemPool masterPool, slavePool, pctgPool;
    
        // load master and slave contigs in memory
        masterPool.loadPool(masterFasta);
        slavePool.loadPool(slaveFasta);
        
        // load reference data (index => contig name)
        BamTools::RefVector mcRef = inBamM.GetReferenceData();
        BamTools::RefVector scRef = inBamS.GetReferenceData();
        
        std::cerr << "[building paired contigs]" << std::endl;
        // build paired contigs (threaded)
        ThreadedBuildPctg tbp( pcblocks, &pctgPool, &masterPool, &slavePool, &mcRef, &scRef );
        std::pair< std::list<PairedContig>, std::vector<bool> > result = tbp.run(threadsNum);
        
        IdType pctgNum((result.first).size());
        std::cerr << "[paired contigs built: " << pctgNum << "]" << std::endl;
        
        slavePool.clear(); // slave contigs pool is no more needed

        std::list<IdType> ctgIds;
        
        for(IdType i = 0; i < result.second.size(); i++) 
            if( !result.second.at(i) ) ctgIds.push_back(i);

        // add unused contigs in paired contigs pool
        generateSingleCtgPctgs( result.first, ctgIds, &pctgPool, &masterPool, &mcRef, pctgNum);

        // save paired contig pool to disk
        std::cerr << "[writing paired contigs on file: " << getPathBaseName( outFilePrefix + ".fasta" ) << "]" << std::endl;
        pctgPool.savePool( outFilePrefix + ".fasta" );
        
        // save paired contigs descriptors to file
        std::cerr << "[writing paired contigs descriptors on file: " << getPathBaseName( outFilePrefix + ".pctgs" ) << "]" << std::endl;
        result.first.sort( orderPctgsByName );
        std::fstream pctgDescFile( (outFilePrefix + ".pctgs").c_str(), std::fstream::out );
        writePctgDescriptors( pctgDescFile, result.first );
        pctgDescFile.close();
        
        inBamM.Close();
        inBamS.Close();
    }
    
    t2 = time(NULL);
    std::cerr << "[execution time: " << formatTime(t2-t1) << "]" << std::endl;
    
    return 0;
}
