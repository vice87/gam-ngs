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
    std::pair< IdType,IdType > ctgid = std::make_pair( block.getMasterFrame().getAssemblyId(), block.getMasterFrame().getContigId() );
    return (this->_masterMap.find(ctgid))->second;
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
typename PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::Vertex 
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::getSlaveVertex( const Block &block ) const
{
    std::pair< IdType,IdType > ctgid = std::make_pair( block.getSlaveFrame().getAssemblyId(), block.getSlaveFrame().getContigId() );
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
        std::pair< IdType,IdType > masterCtgId = std::make_pair( (block->getMasterFrame()).getAssemblyId(), (block->getMasterFrame()).getContigId() );
        std::pair< IdType,IdType > slaveCtgId = std::make_pair( (block->getSlaveFrame()).getAssemblyId(), (block->getSlaveFrame()).getContigId() );
                
        add_edge( this->_masterMap[masterCtgId], this->_slaveMap[slaveCtgId], *this );
    }
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
void 
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::initVertexLabels( const std::vector<Block>& blocks )
{
    // collect master and slave contigs
    std::set< std::pair<IdType,IdType> > masterCtgs;
    std::set< std::pair<IdType,IdType> > slaveCtgs;
    
    //IdType masterCtgId;
    //IdType slaveCtgId;
    std::pair< IdType,IdType > masterCtgId;
    std::pair< IdType,IdType > slaveCtgId;
    
    typename std::vector<Block>::const_iterator block;
    for( block = blocks.begin(); block != blocks.end(); block++ )
    {
        masterCtgId = std::make_pair( (block->getMasterFrame()).getAssemblyId(), (block->getMasterFrame()).getContigId() );
        slaveCtgId = std::make_pair( (block->getSlaveFrame()).getAssemblyId(), (block->getSlaveFrame()).getContigId() );
        
        masterCtgs.insert( masterCtgId );
        slaveCtgs.insert( slaveCtgId );
    }
    
    this->_vertexToCtg.resize( masterCtgs.size() + slaveCtgs.size() );
    
    size_t i = 0;
    
    // insert into the vector the master contig IDs first.
    typename std::set< std::pair<IdType,IdType> >::iterator label;
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


template <class VERTEX_PROP, class EDGE_WEIGHT>
std::ostream&
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::writeGraphviz(std::ostream& os)
{    
    os << "graph AssemblyGraph {" << std::endl;
    //os << "   rankdir=LR;" << std::endl;
    
    VertexIterator vbegin,vend;
    boost::tie(vbegin,vend) = boost::vertices(*this);
    
    for (VertexIterator v=vbegin; v!=vend; v++) 
    {
        os << "\t" << *v << "[label=\"<" << this->_vertexToCtg[*v].first << "," << this->_vertexToCtg[*v].second << ">\"];" << std::endl;
    }
    
    EdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::edges(*this);
    
    for( EdgeIterator e = ebegin; e != eend; e++ )
    {
        os << "\t" << boost::source(*e,*this) << "--" << boost::target(*e,*this) << "[color=black];" << std::endl;
    }
    
    os << "}" << std::endl;
    
    return os;    
}


#endif // PAIRED_GRAPH_CODE_