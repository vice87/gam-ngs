#ifndef PAIRED_GRAPH_CODE_
#define PAIRED_GRAPH_CODE_

#include <vector>
#include <iterator>

#include "graphs/PairedGraph.hpp"

template <class VERTEX_PROP, class EDGE_WEIGHT>
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::PairedGraph() {}


template <class VERTEX_PROP, class EDGE_WEIGHT>
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::PairedGraph(
        const PairedGraph<VERTEX_PROP,EDGE_WEIGHT> &orig)
                : Graph((Graph)orig), _vertexToCtg(orig._vertexToCtg),
                  _masterMap(orig._masterMap), _slaveMap(orig._slaveMap),
                  _firstSlaveVertex(orig._firstSlaveVertex) {}


template <class VERTEX_PROP, class EDGE_WEIGHT>
const PairedGraph<VERTEX_PROP,EDGE_WEIGHT>& 
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::operator =(
        const PairedGraph<VERTEX_PROP,EDGE_WEIGHT> &orig)
{
    *((Graph *)this) = *((Graph *)&orig);
    this->_vertexToCtg = orig._vertexToCtg;
    this->_masterMap = orig._masterMap;
    this->_slaveMap = orig._slaveMap;
    this->_firstSlaveVertex = orig._firstSlaveVertex;
    
    return *this;
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
typename PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::Vertex 
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::getMasterVertex( const Block &block ) const
{
    IdType ctgid = block.getMasterFrame().getContigId();
    return (this->_masterMap.find(ctgid))->second;
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
typename PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::Vertex 
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::getSlaveVertex( const Block &block ) const
{
    IdType ctgid = block.getSlaveFrame().getContigId();
    return (this->_slaveMap.find(ctgid))->second;
}

template <class VERTEX_PROP, class EDGE_WEIGHT>
bool 
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::isMasterNode( const Vertex &node ) const
{
    return (node < this->_firstSlaveVertex);
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
bool 
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::isSlaveNode( const Vertex &node ) const
{
    return (node >= this->_firstSlaveVertex);
}



template <class VERTEX_PROP, class EDGE_WEIGHT>
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::PairedContigGraph() {}


template <class VERTEX_PROP, class EDGE_WEIGHT>
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::PairedContigGraph( 
        const PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT> &orig ) 
                : PairedGraph<VERTEX_PROP,EDGE_WEIGHT>(orig) {}


template <class VERTEX_PROP, class EDGE_WEIGHT>
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::PairedContigGraph( const std::vector<Block> &blocks ) 
{
    this->initGraph( blocks );
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
void 
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::initGraph( const std::vector<Block> &blocks )
{
    this->initVertexLabels(blocks);
    
    // add an edge connected the master and slave contigs vertices of a block.
    typename std::vector<Block>::const_iterator block;
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        IdType masterCtgId = (block->getMasterFrame()).getContigId();
        IdType slaveCtgId = block->getSlaveFrame().getContigId();
                
        add_edge( this->_masterMap[masterCtgId], this->_slaveMap[slaveCtgId], *this );
    }
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
void 
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::initVertexLabels( const std::vector<Block>& blocks )
{
    // collect master and slave contigs
    std::set< IdType > masterCtgs, slaveCtgs;
    
    IdType masterCtgId;
    IdType slaveCtgId;
    
    typename std::vector<Block>::const_iterator block;
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        masterCtgId = (block->getMasterFrame()).getContigId();
        slaveCtgId = (block->getSlaveFrame()).getContigId();
        
        masterCtgs.insert( masterCtgId );
        slaveCtgs.insert( slaveCtgId );
    }
    
    this->_vertexToCtg.resize( masterCtgs.size() + slaveCtgs.size() );
    
    size_t i = 0;
    
    // insert into the vector the master contig IDs first.
    typename std::set< IdType >::iterator label;
    for( label = masterCtgs.begin(); label != masterCtgs.end(); label++ )
    {
        this->_vertexToCtg.at(i) = *label;
        this->_masterMap[*label] = i;
        add_vertex(*this);
        i++;
    }
    
    this->_firstSlaveVertex = i; // index of the first slave contig vertex
    
    // insert into the vector the remaining contig IDs (of slave assembly).
    for( label = slaveCtgs.begin(); label != slaveCtgs.end(); label++ )
    {
        this->_vertexToCtg.at(i) = *label;
        this->_slaveMap[*label] = i;
        add_vertex(*this);
        i++;
    }
}

#endif // PAIRED_GRAPH_CODE_