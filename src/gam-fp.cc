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

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "types.hpp"
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

void printGamCmdError( char** argv )
{
    std::cerr << std::endl << "Usage:   " << getPathBaseName(argv[0]) << " <command> [options]\n";
    std::cerr << "\n";
    std::cerr << "Command: " << "blocks" << "\t" << "generate blocks file\n"
              << "         " << "merge" << "\t" << "merge an assembly with a set of fosmid-pool assemblies\n"
              << std::endl;
}

void printBlocksCmdError( char** argv )
{
    std::cerr << "Usage: " << getPathBaseName(argv[0]) << " blocks "
            << "<master.coordinate-sorted.bam> " 
            << "<fosmid.pools.bams.paths.txt> "
            << "<min-block-size> <output.filename>"
            << std::endl;
}

void printMergeCmdError( char** argv )
{
    std::cerr << "Usage: " << getPathBaseName(argv[0]) << " merge"
            << " [-t threads]"
            << " [-b min-block-size]"
            << " <input.blocks>"
            << " <master.coordinate-sorted.bam>" 
            << " <fosmid.pools.bams.paths.txt>"
            << " <master.fasta>"
            << " <fosmid.pools.fasta.paths.txt>"
            << " <output.filename>"
            << std::endl;
}

std::vector<std::string> loadFileList( std::string file_list )
{
    std::ifstream ifs( file_list.c_str(), std::ifstream::in );
    std::vector< std::string > output;
    
    while (ifs.good())
    {
        std::string file_str;
        getline(ifs,file_str);
        
        if( file_str != "" ) output.push_back( file_str );
    }
    
    return output;
}

/*
 * 
 */
int main(int argc, char** argv) 
{
    std::string action;                 // command
    std::string outFilePrefix;          // output file prefix
    std::string inBlocksFile;           // input blocks file
    
    std::string masterBamFile;          // file BAM (Master)
    std::vector<std::string> fpBamFiles;  // list of fosmid-pool bam files
    
    int minBlockSize = 1;               // dimensione minima (# reads) per i blocchi
    
    std::string masterFasta;            // FASTA file of the master assembly
    std::vector<std::string> fpFastas;    // FASTA files of the fosmid-pool assemblies
    
    size_t threadsNum = 1;              // Numero di thread da usare per il weaving
    
    time_t t1,t2;
    t1 = time(NULL);
    
    int i=1; // parameters index
        
    if( argc < 2 )
    {
        printGamCmdError(argv);       
        return 1;
    }
    
    action = argv[i]; // get command string
    i++;
    
    if(action.compare("blocks") != 0 && action.compare("merge") != 0 )
    {
        printGamCmdError(argv);
        return 1;
    }   
    
    /* Command to build the blocks */
    if( action.compare("blocks") == 0 )
    {
        if( argc != 6 ){ printBlocksCmdError(argv); return 1; }
        
        masterBamFile = argv[2];
        fpBamFiles = loadFileList( std::string(argv[3]) );
        
        minBlockSize = atoi(argv[4]);
        outFilePrefix = argv[5];
        
        if( minBlockSize < 1 ) minBlockSize = 1;
        
        struct stat st;
        
        /* Check for master bam files existence */
        if( stat(masterBamFile.c_str(),&st) != 0 )
        {
            std::cerr << "Master BAM file doesn't exist:\n" << masterBamFile << std::endl;
            return 1;
        }
        
        /* Check for fosmid-pool bam files existence */
        for( std::vector<std::string>::iterator fp_file = fpBamFiles.begin(); fp_file != fpBamFiles.end(); fp_file++ )
        {
            if( stat(fp_file->c_str(),&st) != 0 ) 
            {
                std::cerr << "[error: fosmid-pool bam file " << *fp_file << "doesn't exist]" << std::endl;
                return 1;
            }
        }
        
        if( fpBamFiles.size() == 0 ) // no fosmid-pool bam file exists.
        {
            std::cerr << "[error: no fosmid-pool bam has been given in input]" << std::endl;
            return 1;
        }
        
        BamReader masterBam;
        
        masterBam.Open( masterBamFile );
        if( stat((masterBamFile + ".bai").c_str(),&st) == 0 ) masterBam.OpenIndex( masterBamFile + ".bai" );

        std::cout << "[loading reads in memory]" << std::endl;
        
        /* load only useful reads of the slave assembly */
        sparse_hash_map< std::string, Read > readMap;
        readMap.set_deleted_key("");
        
        // if read-name sorted bams are not provided, load each mapped read.
        Read::loadReadsMap( masterBam, readMap );
        
        std::cout << "[finding blocks]" << std::endl;
        
        int fp_id = 0;
        std::vector<Block> blocks;
        //std::ofstream output( (outFilePrefix + ".ai").c_str() );
        
        for( std::vector<std::string>::const_iterator fp_file = fpBamFiles.begin(); fp_file != fpBamFiles.end(); fp_file++ )
        {
            BamReader fpBam;
            fpBam.Open(*fp_file);
            
            Block::findBlocks( blocks, masterBam, fpBam, minBlockSize, readMap, fp_id );
            //output << getBaseFileName( getPathBaseName( *fp_file ) ) << "\t" << fp_id << std::endl;
            
            fp_id++;
        }
        
        //output.close();
        
        std::cout << "[writing blocks on file: " << getPathBaseName( outFilePrefix + ".blocks" ) << "]" << std::endl;
        Block::writeBlocks( outFilePrefix + ".blocks" , blocks ); // write blocks to file
        
        masterBam.Close();
    }
    
    /* Command to merge assemblies */
    if( action.compare("merge") == 0 )
    {
        while( i < argc && argv[i][0] == '-') // finchè ci sono parametri opzionali da caricare
        {
            std::string opt = argv[i++];
            
            if( opt == "-t" ) // number of threads
            {
                if( i >= argc ){ printMergeCmdError(argv); return 1; }
                threadsNum = atoi(argv[i++]);
                if(threadsNum < 1) threadsNum = 1;
            }
            else if( opt == "-b" )
            {
                if( i >= argc ){ printMergeCmdError(argv); return 1; }
                minBlockSize = atoi(argv[i++]);
            }
            else
            {
                std::cerr << opt << " is not a valid option." << std::endl;
                return 1;
            }
        }
        
        if( argc-i != 6 ) // controllo del numero di argomenti obbligatori del comando "merge"
        {
            printMergeCmdError(argv); 
            return 1;
        }
        
        inBlocksFile = argv[i++];
        
        masterBamFile = argv[i++];
        fpBamFiles = loadFileList( std::string(argv[i++]) );
        
        masterFasta = argv[i++];
        fpFastas = loadFileList( std::string(argv[i++]) );
        
        outFilePrefix = argv[i++];
        
        struct stat st;
        
        /* Check blocks file existence */
        
        if( stat(inBlocksFile.c_str(),&st) != 0 ){ std::cerr << "The file " << inBlocksFile << " doesn't exist." << std::endl; return 1; }
        
        /* Check BAM files existence */
        
        if( stat(masterBamFile.c_str(),&st) != 0 )
        {
            std::cerr << "File " << masterBamFile << " (BAM Master CoordinateSorted) doesn't exist." << std::endl;
            return 1;
        }
        
        /* Check FASTA files existence */
        
        if( stat(masterFasta.c_str(),&st) != 0 )
        {
            std::cerr << "File " << masterFasta << " (Master FASTA) doesn't exist." << std::endl;
            return 1;
        }
        
        std::cout << "[loading blocks]" << std::endl;
        std::vector<Block> blocks = Block::readBlocks( inBlocksFile, minBlockSize );
        std::cout << "[blocks loaded: " << blocks.size() << "]" << std::endl;
        
        // keep only blocks between contigs that share at least minEvidence blocks.
        //std::vector< Block > outBlocks = filterBlocksByPairingEvidences( blocks, minEvidence );
        
        std::cout << "[filtering blocks by frames overlap]" << std::endl;
        // remove adjacent blocks if their frames overlap
        // blocks = Block::filterBlocksByOverlaps( blocks );
        blocks = Block::filterBlocksByCoverage( blocks );
        std::cout << "[blocks filtered: " << blocks.size() << "]" << std::endl;
        
        // create the graph of assemblies, remove cycles and keep remaining blocks.
        std::cout << "[removing cycles from assemblies graph]" << std::endl;
        std::list< std::vector<Block> > pcblocks = partitionBlocks( blocks );
        
        std::cout << "[loading reference data]" << std::endl;
        
        //BamReader masterBam;
        //std::vector< BamReader > fpBams( fpBamFiles.size() );
        //std::cout << "[opening BAM files]" << std::endl;
        
        BamReader bamReader;
        
        // load master bam sequences
        bamReader.Open( masterBamFile );
        BamTools::RefVector mcRef = bamReader.GetReferenceData();
        
        // load fosmid-pools bams sequences
        std::vector< BamTools::RefVector > scRef( fpBamFiles.size() );
        
        for( int i=0; i < fpBamFiles.size(); i++ )
        {
            bamReader.Open( fpBamFiles[i] );
            scRef[i] = bamReader.GetReferenceData();
        }
        
        bamReader.Close();
        
        std::map< std::string, int32_t > masterContigs;
        std::vector< std::map<std::string,int32_t> > slaveCtgsVect( fpBamFiles.size() );
        
        BamTools::RefVector::const_iterator ref_iter;
        
        for( ref_iter = mcRef.begin(); ref_iter != mcRef.end(); ref_iter++ ) 
            masterContigs[ ref_iter->RefName ] = ref_iter->RefLength;
        
        for( int i=0; i < slaveCtgsVect.size(); i++ )
        {
            for( ref_iter = scRef[i].begin(); ref_iter != scRef[i].end(); ref_iter++ )
                slaveCtgsVect[i][ ref_iter->RefName ] = ref_iter->RefLength;
        }
        
        std::cout << "[loading contigs in memory]" << std::endl;
        
        ExtContigMemPool masterPool, slavePool(scRef.size());
        HashContigMemPool pctgPool;
        // load master and slave contigs in memory
        masterPool.loadPool( 0, masterFasta, masterContigs );
        for( int i=0; i < slaveCtgsVect.size(); i++ ) slavePool.loadPool( i, fpFastas[i], slaveCtgsVect[i] );
        
        std::cout << "[building paired contigs]" << std::endl;
        // build paired contigs (threaded)
        ThreadedBuildPctg tbp( pcblocks, &pctgPool, &masterPool, &slavePool, &mcRef, &scRef );
        std::pair< std::list<PairedContig>, std::vector<bool> > result = tbp.run(threadsNum);
        
        IdType pctgNum((result.first).size());
        std::cout << "[paired contigs built: " << pctgNum << "]" << std::endl;
        
        //slavePool.clear(); // slave contigs pool is no more needed
        slavePool.resize(0);

        std::list<IdType> ctgIds;
                
        for(IdType i = 0; i < result.second.size(); i++) 
            if( !result.second.at(i) ) ctgIds.push_back(i);

        // add unused contigs in paired contigs pool
        generateSingleCtgPctgs( result.first, ctgIds, &masterPool, &mcRef, pctgNum);

        // save paired contig pool to disk
        std::cout << "[writing paired contigs on file: " << ( outFilePrefix + ".fasta" ) << "]" << std::endl << std::flush;
        //pctgPool.savePool( outFilePrefix + ".fasta" );
        result.first.sort( orderPctgsByName );
        std::ofstream outFasta((outFilePrefix + ".fasta").c_str(),std::ios::out);
        std::list< PairedContig >::const_iterator i;
        for( i = result.first.begin(); i != result.first.end(); i++ ) outFasta << *i << std::endl;
        outFasta.close();
        
        // save paired contigs descriptors to file
        std::cout << "[writing paired contigs descriptors on file: " << ( outFilePrefix + ".pctgs" ) << "]" << std::endl << std::flush;
        std::fstream pctgDescFile( (outFilePrefix + ".pctgs").c_str(), std::fstream::out );
        writePctgDescriptors( pctgDescFile, result.first, mcRef, scRef );
        pctgDescFile.close();
        
        // save IDs of (slave) contigs NOT merged
        std::cout << "[writing unused slave contigs on file: " << ( outFilePrefix + ".unused" ) << "]" << std::endl << std::flush;
        std::fstream unusedCtgsFile( (outFilePrefix + ".unused").c_str(), std::fstream::out );
        std::vector< std::vector<bool> > usedCtgs;
        for( size_t i=0; i < scRef.size(); i++ ) usedCtgs.push_back( std::vector<bool>( scRef[i].size(), false ) );
        
        for( std::list< PairedContig >::const_iterator pctg = result.first.begin(); pctg != result.first.end(); pctg++ ) 
        {
            typedef std::map< std::pair<IdType,IdType>, ContigInPctgInfo > ContigInfoMap;
            ContigInfoMap slaveCtgs = pctg->getSlaveCtgMap();
            
            for( ContigInfoMap::const_iterator ctg = slaveCtgs.begin(); ctg != slaveCtgs.end(); ctg++ ) 
                usedCtgs[(ctg->first).first][(ctg->first).second] = true;
        }
        
        for( unsigned int i = 0; i < usedCtgs.size(); i++ )
        {
            for( size_t j=0; j < usedCtgs[i].size(); j++ )
                if( !usedCtgs[i][j] ){ unusedCtgsFile << i << "\t" << scRef[i][j].RefName << "\n"; }
        }
        
        unusedCtgsFile.close();
    }
    
    t2 = time(NULL);
    std::cout << "[execution time: " << formatTime(t2-t1) << "]" << std::endl << std::flush;
    
    return 0;
}
