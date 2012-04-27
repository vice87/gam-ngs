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
partitionBlocks( const std::vector< Block > &blocks, const Options &options );


//! Partitions a vector of blocks by paired contigs.
/*!
 * \param blocks a vector of blocks.
 * \return a vector of partitions (vectors of blocks).
 */
std::vector< std::vector< Block > > 
partitionBlocksByPairedContigs( const std::vector< Block > &blocks );

#endif	/* PARTITIONFUNCTIONS_HPP */

