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
    const BamTools::RefVector *_masterRefVector;
    const std::vector< BamTools::RefVector > *_slaveRefVector;
	const Options &_options;

    // output
    std::list< PairedContig > _pctgList;
    std::vector<bool> _removedCtgs;

    // mutex
    pthread_mutex_t _mutexRemoveCtgId;
    pthread_mutex_t _mutex;

    std::pair< std::vector<Block>, bool > extractNextPctg();

    void memorizeNewPctg( std::list< PairedContig > &pctgList );

public:

    ThreadedBuildPctg(
            const std::list< std::vector<Block> > &blocksList,
            HashContigMemPool *pctgPool,
            const ExtContigMemPool *masterPool,
            const ExtContigMemPool *slavePool,
            const BamTools::RefVector *masterRefVector,
            const std::vector< BamTools::RefVector > *slaveRefVector,
            const Options &options
            );


    IdType readPctgNumAndIncrease();

    std::pair< std::list<PairedContig>, std::vector<bool> >
    run();

    friend void*
    buildPctgThread(void *argv);

};


#endif	/* THREADEDBUILDPCTG_HPP */

