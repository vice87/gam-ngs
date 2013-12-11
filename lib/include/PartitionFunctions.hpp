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

#ifndef PARTITIONFUNCTIONS_HPP
#define	PARTITIONFUNCTIONS_HPP

#include <list>

#include "Options.hpp"
#include "assembly/Block.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "graphs/CompactAssemblyGraph.hpp"

using namespace options;

//! Partitions a vector of blocks.
/*!
 * Returns a list of vectors. Each vector contains only blocks which may
 * be merged together.
 *
 * \param blocks a vector of blocks.
 * \param options options of the application.
 * \return a list of block vectors.
 */
std::list< CompactAssemblyGraph* >
partitionBlocks( const std::list<Block> &blocks );


std::vector<double>
computeZScore( MultiBamReader &multiBamReader, const uint64_t &refID, uint32_t start, uint32_t end );

//! Partitions a list of blocks by paired contigs.
/*!
 * \param blocks list of blocks.
 * \return vector of partitions (list of blocks).
 */
std::vector< std::list<Block> >
partitionBlocksByPairedContigs( const std::list< Block > &blocks );

#endif	/* PARTITIONFUNCTIONS_HPP */