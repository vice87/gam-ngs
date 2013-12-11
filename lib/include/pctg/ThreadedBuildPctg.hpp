/*
 *  This file is part of GAM-NGS.
 *  Copyright (c) 2011 by Riccardo Vicedomini <rvicedomini@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Simone Scalabrin <scalabrin@appliedgenomics.org>,
 *  Lars Arverstad <lars.arvestad@scilifelab.se>,
 *  Alberto Policriti <policriti@appliedgenomics.org>,
 *  Alberto Casagrande <casagrande@appliedgenomics.org>
 *
 *  GAM-NGS is an evolution of a previous work (GAM) done by Alberto Casagrande,
 *  Cristian Del Fabbro, Simone Scalabrin, and Alberto Policriti.
 *  In particular, GAM-NGS has been adapted to work on NGS data sets and it has
 *  been written using GAM's software as starting point. Thus, it shares part of
 *  GAM's source code.
 *
 *  GAM-NGS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GAM-NGS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GAM-NGS.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef THREADEDBUILDPCTG_HPP
#define	THREADEDBUILDPCTG_HPP

#include <pthread.h>
#include <list>

#include "bam/MultiBamReader.hpp"
#include "graphs/CompactAssemblyGraph.hpp"
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

    std::vector< CompactAssemblyGraph* > _graphs;
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
    CompactAssemblyGraph* extractNextPctg();
	IdType readPctgNumAndIncrease();
	void incProcBlocks( uint64_t num, uint64_t tid );

public:

    ThreadedBuildPctg(
            const std::list< CompactAssemblyGraph* > &graphsList,
            const RefSequence &masterRef,
            const RefSequence &slaveRef
	);

    std::list<PairedContig>* run();

	double computeZScore( MultiBamReader &multiBamReader, int32_t ctgId, uint32_t start, uint32_t end, bool isMaster );

    friend void* buildPctgThread(void *argv);
};


#endif	/* THREADEDBUILDPCTG_HPP */