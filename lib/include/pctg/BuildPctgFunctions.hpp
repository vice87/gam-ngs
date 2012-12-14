/*
 * File:   BuildPctgFunctions.hpp
 * Author: vice
 *
 * Created on 4 giugno 2011, 23.29
 */

#ifndef BUILDPCTGFUNCTIONS_HPP
#define	BUILDPCTGFUNCTIONS_HPP

#include <iostream>
#include <vector>
#include <list>
#include <map>

#include <boost/config.hpp>
#include <boost/graph/topological_sort.hpp>

#include "api/BamAux.h"
#include "api/BamReader.h"

#include "graphs/AssemblyGraph.hpp"
#include "graphs/CompactAssemblyGraph.hpp"
#include "pctg/PairedContig.hpp"
#include "pctg/PctgBuilder.hpp"
#include "pctg/ThreadedBuildPctg.hpp"
#include "pool/HashContigMemPool.hpp"
#include "assembly/io_contig.hpp"


std::list< PairedContig >& buildPctg(
        CompactAssemblyGraph &ag,
        PctgBuilder &builder,
        std::list< PairedContig > &pctgList
);


std::list< PairedContig >& buildPctg(
		ThreadedBuildPctg *tbp,
        CompactAssemblyGraph &ag,
        const RefSequence &masterRef,
		const RefSequence &slaveRef,
        std::list< PairedContig > &pctgList
);


void generateSingleCtgPctgs(
        std::list<PairedContig> &pctgList,
        const std::list<int32_t> &ctgIds,
        const RefSequence &masterRef,
        uint64_t &pctgId
);

#endif	/* BUILDPCTGFUNCTIONS_HPP */

