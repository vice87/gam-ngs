/*! 
 * \file OrderingFunctions.hpp
 * \brief Definition of utility functions for blocks ordering.
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
        return (
                 a.getMasterFrame().getAssemblyId() < b.getMasterFrame().getAssemblyId() ||
                (a.getMasterFrame().getAssemblyId() == b.getMasterFrame().getAssemblyId() && a.getMasterFrame().getContigId() < b.getMasterFrame().getContigId()) || 
                (a.getMasterFrame().getAssemblyId() == b.getMasterFrame().getAssemblyId() && a.getMasterFrame().getContigId() == b.getMasterFrame().getContigId() && a.getMasterFrame().getBegin() < b.getMasterFrame().getBegin()) ||
                (a.getMasterFrame().getAssemblyId() == b.getMasterFrame().getAssemblyId() && a.getMasterFrame().getContigId() == b.getMasterFrame().getContigId() && a.getMasterFrame().getBegin() == b.getMasterFrame().getBegin() && a.getMasterFrame().getLength() > b.getMasterFrame().getLength())
               );
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
        return (
                 a.getSlaveFrame().getAssemblyId() < b.getSlaveFrame().getAssemblyId() ||
                (a.getSlaveFrame().getAssemblyId() == b.getSlaveFrame().getAssemblyId() && a.getSlaveFrame().getContigId() < b.getSlaveFrame().getContigId()) || 
                (a.getSlaveFrame().getAssemblyId() == b.getSlaveFrame().getAssemblyId() && a.getSlaveFrame().getContigId() == b.getSlaveFrame().getContigId() && a.getSlaveFrame().getBegin() < b.getSlaveFrame().getBegin()) ||
                (a.getSlaveFrame().getAssemblyId() == b.getSlaveFrame().getAssemblyId() && a.getSlaveFrame().getContigId() == b.getSlaveFrame().getContigId() && a.getSlaveFrame().getBegin() == b.getSlaveFrame().getBegin() && a.getSlaveFrame().getLength() > b.getSlaveFrame().getLength()) 
               );
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
 * \param blocks vector of blocks
 * \param comparison struct for Block objects
 * \return a pair consisting of ordered indices and back indices
 */
template< class BlocksOrderer >
inline std::pair< std::vector<UIntType>, std::vector<UIntType> > 
getOrderedIndices( const std::vector<Block> &blocks, BlocksOrderer orderer )
{
    std::map< Block, UIntType > blockMap;
    
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
inline getOrderedMasterIndices( const std::vector<Block> &blocks )
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
inline getOrderedSlaveIndices( const std::vector<Block> &blocks )
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

