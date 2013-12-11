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

/*!
 * \file PctgBuilder.hpp
 * \brief Definition of PctgBuilder class.
 * \details This file contains the definition of the class that implements the
 * construction of a paired contig.
 */

#ifndef PCTGBUILDER_HPP
#define	PCTGBUILDER_HPP

#include <iostream>
#include <stdexcept>
#include <pthread.h>

#include "api/BamAux.h"
#include "api/BamMultiReader.h"

#include "types.hpp"
#include "alignment/my_alignment.hpp"
#include "assembly/Block.hpp"
#include "assembly/RefSequence.hpp"
#include "pctg/BestCtgAlignment.hpp"
#include "pctg/ContigInPctgInfo.hpp"
#include "pctg/CtgInPctgInfo.hpp"
#include "pctg/PairedContig.hpp"
#include "pctg/constraints_disattended.hpp"
#include "pctg/ThreadedBuildPctg.hpp"
#include "PartitionFunctions.hpp"

#include "graphs/CompactAssemblyGraph.hpp"
#include "pctg/MergeDescriptor.hpp"

#ifndef MIN_HOMOLOGY
#define MIN_HOMOLOGY 95.0
#endif

#ifndef MIN_HOMOLOGY_2
#define MIN_HOMOLOGY_2 90.0
#endif

#ifndef MIN_ALIGNMENT_LEN
#define MIN_ALIGNMENT_LEN 100
#endif

#define DEFAULT_MAX_GAPS 300
#define DEFAULT_MAX_SEARCHED_ALIGNMENT 400000

#ifndef MIN_ALIGNMENT_QUOTIENT
#define MIN_ALIGNMENT_QUOTIENT 0.001
#endif


//! Class implementing a builder of paired contigs.
class PctgBuilder
{

private:
	ThreadedBuildPctg *_tbp;

	const RefSequence *_masterRef;        //!< Reference to the id->name vector of the master contigs.
	const RefSequence *_slaveRef;         //!< Reference to the id->name vector of the slave contigs.

    UIntType _maxAlignment;                             //!< Maximum alignment size
    UIntType _maxPctgGap;                               //!< Maximum paired contig gaps
    UIntType _maxCtgGap;                                //!< Maximum contig gaps

public:
    //! A constructor.
    /*!
     * Creates a paired contig builder.
     * \param masterPool reference to the master contigs pool
     * \param slavePool reference to the slave contigs pool
     * \param masterRefVector reference to the vector id->name of master contigs
     * \param slaveRefVector reference to the vector id->name of slave contigs
     */
    PctgBuilder(
			ThreadedBuildPctg *tbp = NULL,
            const RefSequence *masterRef = NULL,
            const RefSequence *slaveRef = NULL);

    //! Gets a master contig, given its ID.
    /*!
     * \param ctgId pair consisting in master sequence's assembly and contig identifiers.
     * \return reference to the master contig with ID \c ctgId.
     */
    inline const Contig& loadMasterContig(const int32_t ctgId) const
	{
		return *(this->_masterRef->at(ctgId)).Sequence;
	}

    //! Gets a slave contig, given its ID.
    /*!
     * \param ctgId pair consisting in slave sequence's assembly and contig identifiers.
     * \return reference to the slave contig with ID \c ctgId.
     */
    inline const Contig& loadSlaveContig(const int32_t ctgId) const
	{
		return *(this->_slaveRef->at(ctgId)).Sequence;
	}

	//! Adds the first contig to a paired contig.
    /*!
     * The contig must be one of the master assembly.
     * \param pctg a paired contig.
     * \param ctgId a pair consisting in assembly and contig identifiers of the sequence.
     * \return a paired contig consisting of the single (master) contig \c ctgId.
     *
     * \throws std::invalid_argument if the \c pctg is not empty.
     */
    PairedContig& addFirstContigTo(PairedContig &pctg, const int32_t ctgId) const;


    //! Builds a paired contig with a single (master) contig.
    /*!
     * \param pctgId the paired contig's identifier
     * \param ctgId  a (master) contig identifier
     */
    PairedContig initByContig(const IdType &pctgId, const int32_t ctgId) const;

	void appendMasterToPctg( PairedContig &pctg, int32_t id, Contig &ctg, int32_t start, int32_t end, bool rev );
	void appendSlaveToPctg( PairedContig &pctg, int32_t id, Contig &ctg, int32_t start, int32_t end, bool rev );
	void appendBlocksRegionToPctg( PairedContig &pctg, int32_t m_id, Contig &m_ctg, int32_t m_start, int32_t m_end, bool m_rev,
								   int32_t s_id, Contig &s_ctg, int32_t s_start, int32_t s_end, bool s_rev );

	void buildPctgs( std::list<PairedContig> &pctgList, MergeBlockLists &mergeLists );
	void buildPctgs( std::list<PairedContig> &pctgList, std::list<MergeBlock> &ml );

	void splitMergeBlocksByInclusions( MergeBlockLists &ml_in );
	void sortMergeBlocksByDirection( MergeBlockLists &ml );
	void splitMergeBlocksByDirection( MergeBlockLists &ml_in );
	void splitMergeBlocksByAlign( MergeBlockLists &ml_in );
	void alignMergeBlock( const CompactAssemblyGraph &graph, MergeBlock &mb ) const;
    bool getMergePaths( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::vector<MergeBlock> &mbv,
                        MergeBlockLists &merge_paths ) const;
    bool solveForks( CompactAssemblyGraph &graph, std::vector<MergeBlock> &mbv ) const;
    void getNextContigs( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::set<int32_t> &masterCtgs, std::set<int32_t> &slaveCtgs ) const;
    void getPrevContigs( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::set<int32_t> &masterCtgs, std::set<int32_t> &slaveCtgs ) const;


    //! Computes the best alignment between a paired contig and a contig which may be merged.
    /*!
     * \param pctg a paired contig
     * \param pctgInfo
     * \param startPos position in \c pctg from which to find the alignment (where the first block starts)
	 * \param endPos position in \c pctg where the last block ends
     * \param ctg a contig
	 * \param mergeMasterCtg whether \c ctg is a master or a slave contig
     * \param blocks_list list of blocks between the contigs to be merged
     */
    void findBestAlignment(
        BestCtgAlignment &bestAlign,
        Contig &masterCtg,
		uint64_t masterStart,
		uint64_t masterEnd,
		Contig &slaveCtg,
        uint64_t slaveStart,
        uint64_t slaveEnd,
        const std::list<Block>& blocks_list ) const;

	void alignBlocks(
		const Contig &masterCtg,
		const uint64_t &masterStart,
		const Contig &slaveCtg,
		const uint64_t &slaveStart,
		const std::list<Block> &blocks_list,
		std::vector< MyAlignment > &alignments ) const;

	bool is_good( const std::vector<MyAlignment> &align, uint64_t min_align_len = MIN_ALIGNMENT_LEN ) const;
	bool is_good( const MyAlignment &align, uint64_t min_align_len = MIN_ALIGNMENT_LEN ) const;
};

#endif	/* PCTGBUILDER_HPP */

