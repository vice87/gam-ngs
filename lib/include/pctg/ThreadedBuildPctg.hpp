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

#include "pctg/PairedContig.hpp"

void * buildPctgThread(void *argv);

class ThreadedBuildPctg
{

private:

    typedef struct thread_arg
    {
            ThreadedBuildPctg *tbp;
            uint64_t tid;
            std::list< PairedContig > *output;

    } thread_arg_t;

    // input
    const RefSequence& _masterRef;
    const RefSequence& _slaveRef;

    IdType _pctgNum;
    UIntType _pctgsDone;
    UIntType _lastPerc;

    std::vector< std::list<Block>* > _blocksVect;
    uint64_t _nextPctg;

    uint64_t _procBlocks;
    uint64_t _totBlocks;

    // output
    std::list< PairedContig > _pctgList;

    // mutex
    pthread_mutex_t _mutexRemoveCtgId;
    pthread_mutex_t _mutex;
    pthread_mutex_t _mutexProcBlocks;
    pthread_mutex_t _mutexPctgNumInc;

    pthread_mutex_t _mutexMasterBam; // mutex per accedere al BAM master
    pthread_mutex_t _mutexSlaveBam; // mutex per accedere al BAM slave

	// private methods
    std::list<Block>* extractNextPctg();
	IdType readPctgNumAndIncrease();
	void incProcBlocks( uint64_t num, uint64_t tid );

public:

    ThreadedBuildPctg(
            const std::list< std::vector<Block> > &blocksList,
            const RefSequence &masterRef,
            const RefSequence &slaveRef
	);

    std::list<PairedContig>* run();

	double computeZScore( MultiBamReader &multiBamReader, int32_t ctgId, uint32_t start, uint32_t end, bool isMaster );

    friend void* buildPctgThread(void *argv);
};


#endif	/* THREADEDBUILDPCTG_HPP */