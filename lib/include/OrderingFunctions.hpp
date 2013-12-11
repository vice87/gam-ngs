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

#ifndef ORDERINGFUNCTIONS_HPP
#define	ORDERINGFUNCTIONS_HPP

#include <map>
#include <vector>
#include <algorithm>

#include "assembly/Block.hpp"

//! Comparison struct for blocks according to their master frame
struct MasterBlocksOrderer
{
    //! Returns whether a block comes before another one in their master frame
    /*!
     * \param a a block
     * \param b a block
     * \return \c true if a comes before b in their master contig
     */
    bool operator()( const Block &a, const Block &b )
    {
		const Frame& af = a.getMasterFrame();
		const Frame& bf = b.getMasterFrame();

		int32_t a_id = af.getContigId();
		int32_t b_id = bf.getContigId();

		if( a_id < b_id ) return true;
		if( a_id > b_id ) return false;

		int32_t a_beg = af.getBegin();
		int32_t b_beg = bf.getBegin();

		if( a_beg < b_beg ) return true;
		if( a_beg > b_beg ) return false;

		int32_t a_len = af.getLength();
		int32_t b_len = bf.getLength();

		return ( a_len > b_len );
    }
};

//! Comparison struct for blocks according to their slave frame
struct SlaveBlocksOrderer
{
    //! Returns whether a block comes before another one in their slave frame
    /*!
     * \param a a block
     * \param b a block
     * \return \c true if a comes before b in their slave contig
     */
    bool operator()( const Block &a, const Block &b )
    {
		const Frame& af = a.getSlaveFrame();
		const Frame& bf = b.getSlaveFrame();

		int32_t a_id = af.getContigId();
		int32_t b_id = bf.getContigId();

		if( a_id < b_id ) return true;
		if( a_id > b_id ) return false;

		int32_t a_beg = af.getBegin();
		int32_t b_beg = bf.getBegin();

		if( a_beg < b_beg ) return true;
		if( a_beg > b_beg ) return false;

		int32_t a_len = af.getLength();
		int32_t b_len = bf.getLength();

		return ( a_len > b_len );
    }
};


//! Comparison struct for blocks according to their number of reads
struct ReadNumBlocksOrderer
{
    //! Returns whether a block comes before another one in their slave frame
    /*!
     * \param a a block
     * \param b a block
     * \return \c true if a has less reads than b
     */
    bool operator()( const Block &a, const Block &b )
    {
        return a.getReadsNumber() < b.getReadsNumber();
    }
};


//! Returns a pair consisting of the ordered indices and the back indices of a vector of blocks
/*!
 * <ul>
 *      <li> <b>Ordered indices vector</b>: the index in position \c i indicates the index of the <code>i</code>-th block in \c blocks. </li>
 *      <li> <b>Back indices vector</b>: given an index \c i of \c blocks, returns the index of a block in the sorted vector. </li>
 * </ul>
 *
 * \param blocksList list of blocks
 * \param comparison struct for Block objects
 * \return a pair consisting of ordered indices and back indices
 */
template< class BlocksOrderer >
inline std::pair< std::vector<UIntType>, std::vector<UIntType> >
getOrderedIndices( const std::list<Block> &blocksList, BlocksOrderer orderer )
{
    std::map< Block, UIntType > blockMap;
	std::vector< Block > blocks;

	blocks.reserve( blocksList.size() );
	for( std::list<Block>::const_iterator b = blocksList.begin(); b != blocksList.end(); ++b ) blocks.push_back(*b);

    for( UIntType idx = 0; idx < blocks.size(); idx++ ) blockMap[ blocks[idx] ] = idx;

    std::vector< Block > orderedBlocks(blocks);
    std::sort( orderedBlocks.begin(), orderedBlocks.end(), orderer );

    std::vector< UIntType > orderedIdx( blocks.size() );
    std::vector< UIntType > reverseIdx( blocks.size() );

    // fill index vectors
    for( UIntType i = 0; i < blocks.size(); i++ )
    {
        orderedIdx[i] = blockMap[ orderedBlocks[i] ];
        reverseIdx[ orderedIdx[i] ] = i;
    }

    return std::pair< std::vector<UIntType>, std::vector<UIntType> >( orderedIdx, reverseIdx );
}


//! Returns ordered indices and back indices of a vector of blocks according
//! to their master contig
/*!
 * \param blocks vector of blocks
 * \return ordered and back indices
 */
std::pair< std::vector<UIntType>, std::vector<UIntType> >
inline getOrderedMasterIndices( const std::list<Block> &blocks )
{
    MasterBlocksOrderer mbo;
    return getOrderedIndices( blocks, mbo );
}

//! Returns ordered indices and back indices of a vector of blocks according
//! to their slave contig
/*!
 * \param blocks vector of blocks
 * \return ordered and back indices
 */
std::pair< std::vector<UIntType>, std::vector<UIntType> >
inline getOrderedSlaveIndices( const std::list<Block> &blocks )
{
    SlaveBlocksOrderer sbo;
    return getOrderedIndices( blocks, sbo );
}


//! Compare two strings according to samtools sort order
/*!
 * \param a first string
 * \param b second string
 * \return 0 if a == b, -1 if a < b, 1 otherwise
 */
int
inline strnum_cmp(const char *a, const char *b)
{
    char *pa, *pb;
    pa = (char*)a; pb = (char*)b;

    while (*pa && *pb)
    {
        if (isdigit(*pa) && isdigit(*pb))
        {
            long ai, bi;
            ai = strtol(pa, &pa, 10);
            bi = strtol(pb, &pb, 10);

            if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
        }
        else
        {
            if (*pa != *pb) break;
            ++pa; ++pb;
        }
    }

    if (*pa == *pb) return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;

    return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

#endif	/* ORDERINGFUNCTIONS_HPP */

