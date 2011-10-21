/* 
 * File:   ThreadedBuildPctg.hpp
 * Author: vice
 *
 * Created on 12 giugno 2011, 22.35
 */

#ifndef THREADEDBUILDPCTG_HPP
#define	THREADEDBUILDPCTG_HPP

#include <pthread.h>

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
    ContigMemPool *_pctgPool;
    const ContigMemPool *_masterPool;
    const ContigMemPool *_slavePool;
    const BamTools::RefVector *_masterRefVector;
    const BamTools::RefVector *_slaveRefVector;
    
    // output
    std::list< PairedContig > _pctgList;
    std::vector<bool> _removedCtgs;
    
    // mutex
    pthread_mutex_t _mutexRemoveCtgId;
    pthread_mutex_t _mutex;
    
    std::pair< std::vector<Block>, bool > extractNextPctg();
    
    void memorizeNewPctg(const PairedContig &pctg);
    
public:
    
    ThreadedBuildPctg( 
            const std::list< std::vector<Block> > &blocksList,
            ContigMemPool *pctgPool,
            const ContigMemPool *masterPool,
            const ContigMemPool *slavePool,
            const BamTools::RefVector *masterRefVector,
            const BamTools::RefVector *slaveRefVector
            );
    
    
    IdType readPctgNumAndIncrease();
    
    std::pair< std::list<PairedContig>, std::vector<bool> >
    run(const size_t &threadsNum);
    
    friend void*
    buildPctgThread(void *argv);
    
};


#endif	/* THREADEDBUILDPCTG_HPP */

