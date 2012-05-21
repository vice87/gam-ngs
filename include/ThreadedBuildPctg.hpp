/*
 * File:   ThreadedBuildPctg.hpp
 * Author: vice
 *
 * Created on 12 giugno 2011, 22.35
 */

#ifndef THREADEDBUILDPCTG_HPP
#define	THREADEDBUILDPCTG_HPP

#include <pthread.h>
#include <list>

#include <bam/MultiBamReader.hpp>

#include "Options.hpp"
using namespace options;

void *
buildPctgThread(void *argv);


class ThreadedBuildPctg
{

private:

    // input
    IdType _pctgNum;
    UIntType _pctgsDone;
    UIntType _lastPerc;
    std::list< std::vector<Block> > _blocksList;

    HashContigMemPool *_pctgPool;
    const ExtContigMemPool *_masterPool;
    const ExtContigMemPool *_slavePool;

	MultiBamReader& _masterBam;
	std::vector< MultiBamReader >& _slaveBams;

    const BamTools::RefVector *_masterRefVector;
    const std::vector< BamTools::RefVector > *_slaveRefVector;

	const Options &_options;

    // output
    std::list< PairedContig > _pctgList;
    std::vector<bool> _removedCtgs;

    // mutex
    pthread_mutex_t _mutexRemoveCtgId;
    pthread_mutex_t _mutex;

	pthread_mutex_t _mutexMasterBam; // mutex per accedere al BAM master
	std::vector< pthread_mutex_t > _mutexSlaveBams; // mutex per accedere ai BAM slave

    std::pair< std::vector<Block>, bool > extractNextPctg();

    void memorizeNewPctg( std::list< PairedContig > &pctgList );

public:

    ThreadedBuildPctg(
            const std::list< std::vector<Block> > &blocksList,
            HashContigMemPool *pctgPool,
            const ExtContigMemPool *masterPool,
            const ExtContigMemPool *slavePool,
			MultiBamReader &masterBam,
			std::vector< MultiBamReader > &slaveBams,
            const BamTools::RefVector *masterRefVector,
            const std::vector< BamTools::RefVector > *slaveRefVector,
            const Options &options
            );


    IdType readPctgNumAndIncrease();

    std::pair< std::list<PairedContig>, std::vector<bool> >  run();

	double computeZScore( const std::pair<uint64_t,uint64_t> &ctgId, uint32_t start, uint32_t end, bool isMaster );
	double computeReadCoverage( std::pair<uint64_t,uint64_t> &ctgId, uint32_t start, uint32_t end, bool isMaster );

    friend void* buildPctgThread(void *argv);
};


#endif	/* THREADEDBUILDPCTG_HPP */

