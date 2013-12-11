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
 * \file PairingEvidencesGraph.hpp
 * \brief Definition of PairingEvidencesGraph class.
 * \details This class construct a PairedContigGraph whose edges are weighted
 * according to the number of blocks on two contigs.
 */

#ifndef PAIRINGEVIDENCESGRAPH_HPP
#define	PAIRINGEVIDENCESGRAPH_HPP

#include "graphs/PairedGraph.code.hpp"

//! PairingEvidencesGraph class.
/*!
 * Extends the PairedContigGraph class, adding a weight to each edge.
 */
class PairingEvidencesGraph : public PairedContigGraph<>
{

public:
    typedef PairedContigGraph<>::Edge Edge; //!< edge descriptor type
    typedef PairedContigGraph<>::Vertex Vertex; //!< vertex descriptor type

private:
    //! Adds edge weights.
    /*!
     * For each block, the weight of the edge connecting its master and slave contigs
     * is increased by 1.
     * \param blocks a list of blocks.
     */
    void addEdgeWeights( const std::list<Block>& blocks );

public:
    //! A constructor.
    /*!
     * \param blocks a list of blocks.
     */
    PairingEvidencesGraph( const std::list<Block> &blocks );

};

//! Discards blocks whose master and slave contigs have a pairing evidence below a certain threshold.
/*!
 * \param blocks a vector of blocks.
 * \param minPairEvid minimum number of blocks between two contigs.
 * \return a filtered vector of blocks.
 */
//std::vector<Block> filterBlocksByPairingEvidences( const std::vector<Block> &blocks, const int minPairEvid = 1 );

void getSingleLinkBlocks(
	const std::list<Block> &blocks,
	std::set< std::pair<int32_t,int32_t> > &slb
);

#endif	/* PAIRINGEVIDENCESGRAPH_HPP */

