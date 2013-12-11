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

