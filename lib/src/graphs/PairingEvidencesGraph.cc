
#include "graphs/PairingEvidencesGraph.hpp"

PairingEvidencesGraph::PairingEvidencesGraph(const std::vector<Block> &blocks) 
        : PairedContigGraph< >(blocks)
{
    this->addEdgeWeights(blocks);
}

void 
PairingEvidencesGraph::addEdgeWeights( const std::vector<Block>& blocks )
{
    std::pair<IdType,IdType> masterCtgId, slaveCtgId;
    UIntType weight=0;
    
    std::vector<Block>::const_iterator block;
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        masterCtgId.first = (block->getMasterFrame()).getAssemblyId();
        masterCtgId.second = (block->getMasterFrame()).getContigId();
        slaveCtgId.first = (block->getSlaveFrame()).getAssemblyId();
        slaveCtgId.second = (block->getSlaveFrame()).getContigId();
        
        PairingEvidencesGraph::Edge e = 
                boost::edge(this->_masterMap[masterCtgId], this->_slaveMap[slaveCtgId], *this ).first;
        
        put( boost::edge_weight_t(), *this, e, weight );
    }
    
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        masterCtgId.first = (block->getMasterFrame()).getAssemblyId();
        masterCtgId.second = (block->getMasterFrame()).getContigId();
        slaveCtgId.first = (block->getSlaveFrame()).getAssemblyId();
        slaveCtgId.second = (block->getSlaveFrame()).getContigId();
        
        PairingEvidencesGraph::Edge e = 
                boost::edge(this->_masterMap[masterCtgId], this->_slaveMap[slaveCtgId], *this ).first;
        
        weight = get( boost::edge_weight_t(), *this, e );
        weight++;
        put( boost::edge_weight_t(), *this, e, weight);
    }
}


std::vector<Block> filterBlocksByPairingEvidences( const std::vector<Block> &blocks, const int minPairEvid )
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
}