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

#ifndef MERGEDESCRIPTOR_HPP_
#define MERGEDESCRIPTOR_HPP_

#include "pctg/BestCtgAlignment.hpp"
#include "assembly/contig.hpp"

typedef enum { LINEAR_MERGE, FORK_MERGE } __attribute__((packed)) MergeType;
typedef enum { MIS_MASTER, MIS_SLAVE, REPEAT, UNKNOWN } __attribute__((packed)) ForkType;

struct MergeBlock
{
	uint64_t vertex;

	int32_t m_id;
	int32_t m_start;
	int32_t m_end;

	int32_t s_id;
	int32_t s_start;
	int32_t s_end;

    bool valid;

	bool align_rev;
	bool align_ok;

	// relative to orginal strand of contigs
	bool m_ltail;
	bool m_rtail;
	bool s_ltail;
	bool s_rtail;

	// relative to order of merge
	bool ext_slave_next;
	bool ext_slave_prev;

	bool m_rev;
	bool s_rev;
};

struct MergeDescriptor
{
    MergeType mergeType;
    ForkType forkType;

	int32_t masterId;
	int32_t slaveId;

	std::pair<int64_t,int64_t> masterMainAlign;
	std::pair<int64_t,int64_t> slaveMainAlign;

    std::pair<int64_t,int64_t> masterPos;
	std::pair<int64_t,int64_t> slavePos;

	std::pair<int64_t,int64_t> masterBlocks;
	std::pair<int64_t,int64_t> slaveBlocks;

	BestCtgAlignment *align;
};

typedef std::list< std::list<MergeDescriptor> > MergeDescriptorLists;
typedef std::list< std::list<MergeBlock> > MergeBlockLists;


struct MergeStruct
{
	int32_t id;
	Contig *ctg;
	bool reversed;
	int64_t pos;

	std::pair<int64_t,int64_t> limits;
};

#endif //MERGEDESCRIPTOR_HPP_
