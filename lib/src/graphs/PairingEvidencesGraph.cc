
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