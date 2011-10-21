/* 
 * File:   BuildPctgFunctions.hpp
 * Author: vice
 *
 * Created on 4 giugno 2011, 23.29
 */

#ifndef BUILDPCTGFUNCTIONS_HPP
#define	BUILDPCTGFUNCTIONS_HPP

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
#include "pool/ContigMemPool.hpp"
#include "assembly/io_contig.hpp"


std::list< PairedContig > buildPctg( 
        IdType &pctgId,
        ThreadedBuildPctg *tbp,
        const AssemblyGraph &ag, 
        const PctgBuilder &builder );


std::list< PairedContig > buildPctg( 
        const AssemblyGraph &ag,
        const ContigMemPool *masterPool, 
        const ContigMemPool *slavePool, 
        const BamTools::RefVector *masterRefVector,
        const BamTools::RefVector *slaveRefVector,
        IdType &pctgId,
        ThreadedBuildPctg *tbp );


void generateSingleCtgPctgs(
        std::list<PairedContig> &pctgList,
        const std::list<IdType> &ctgIds,
        ContigMemPool *pctgPool,
        ContigMemPool *masterPool,
        BamTools::RefVector *masterRefVector,
        IdType &pctgId
);

#endif	/* BUILDPCTGFUNCTIONS_HPP */

