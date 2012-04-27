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

#include "ThreadedBuildPctg.hpp"

#include "graphs/AssemblyGraph.hpp"
#include "pctg/PairedContig.hpp"
#include "pctg/PctgBuilder.hpp"
#include "pool/HashContigMemPool.hpp"
#include "assembly/io_contig.hpp"


std::list< PairedContig >& buildPctg(
        const AssemblyGraph &ag,
        const PctgBuilder &builder,
        std::list< PairedContig > &pctgList,
		const Options &options );


std::list< PairedContig >& buildPctg(
        const AssemblyGraph &ag,
        const ExtContigMemPool *masterPool,
        const ExtContigMemPool *slavePool,
        const BamTools::RefVector *masterRefVector,
        const std::vector<BamTools::RefVector> *slaveRefVector,
        std::list< PairedContig > &pctgList,
		const Options &options );


void generateSingleCtgPctgs(
        std::list<PairedContig> &pctgList,
        const std::list<IdType> &ctgIds,
        //HashContigMemPool *pctgPool,
        ExtContigMemPool *masterPool,
        BamTools::RefVector *masterRefVector,
        IdType &pctgId
);

#endif	/* BUILDPCTGFUNCTIONS_HPP */

