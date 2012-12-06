/*!
 * \file PartitionFunctions.hpp
 * \brief Definition of utility functions for partitioning blocks which may be
 * merged together.
 */

#ifndef PARTITIONFUNCTIONS_HPP
#define	PARTITIONFUNCTIONS_HPP

#include <list>
#include "assembly/Block.hpp"
#include "Options.hpp"

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
std::list< std::vector< Block > >
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