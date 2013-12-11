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

#include "graphs/PairingEvidencesGraph.hpp"

PairingEvidencesGraph::PairingEvidencesGraph(const std::list<Block> &blocks)
        : PairedContigGraph< >(blocks)
{
    this->addEdgeWeights(blocks);
}

void
PairingEvidencesGraph::addEdgeWeights( const std::list<Block> &blocks )
{
    int32_t masterCtgId, slaveCtgId;
    uint64_t weight=0;

	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
        masterCtgId = b->getMasterId();
        slaveCtgId = b->getSlaveId();

        PairingEvidencesGraph::Edge e = boost::edge( _masterMap[masterCtgId], _slaveMap[slaveCtgId], *this ).first;

        put( boost::edge_weight_t(), *this, e, weight );
    }

    for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
		masterCtgId = b->getMasterId();
		slaveCtgId = b->getSlaveId();

        PairingEvidencesGraph::Edge e = boost::edge( _masterMap[masterCtgId], _slaveMap[slaveCtgId], *this ).first;

        weight = get( boost::edge_weight_t(), *this, e );
        weight++;

        put( boost::edge_weight_t(), *this, e, weight); // update edge weight
    }
}


/*std::vector<Block> filterBlocksByPairingEvidences( const std::vector<Block> &blocks, const int minPairEvid )
{
    PairingEvidencesGraph peg(blocks);

    UIntType outBlockNum = 0;

    // Count blocks with pairing evidences >= minPairEvid
    std::vector<Block>::const_iterator block;
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        PairingEvidencesGraph::Edge e = boost::edge(peg.getMasterVertex(*block), peg.getSlaveVertex(*block), peg).first;

        if( get(boost::edge_weight_t(), peg, e) >= minPairEvid ) outBlockNum++;
    }

    // Fill the output vector with valid blocks
    std::vector<Block> outBlocks( outBlockNum );
    outBlockNum = 0;
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        PairingEvidencesGraph::Edge e = boost::edge(peg.getMasterVertex(*block), peg.getSlaveVertex(*block), peg).first;

        if( get(boost::edge_weight_t(), peg, e) >= minPairEvid ) outBlocks[ outBlockNum++ ] = *block;
    }

    return outBlocks;
}*/


void getSingleLinkBlocks(
	const std::list<Block> &blocks,
	std::set< std::pair<int32_t,int32_t> > &slb )
{
	PairingEvidencesGraph peg(blocks);

	int32_t mid, sid;
	PairingEvidencesGraph::Vertex mv, sv;

	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
	{
		mv = peg.getMasterVertex(*b);
		sv = peg.getSlaveVertex(*b);

		if( boost::out_degree(mv,peg) == 1 || boost::out_degree(sv,peg) == 1 )
			slb.insert( std::make_pair( b->getMasterId(), b->getSlaveId() ) );
	}
}