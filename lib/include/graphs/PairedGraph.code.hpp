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
    int32_t ctgid = block.getMasterId();
    return (this->_masterMap.find(ctgid))->second;
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
typename PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::Vertex
PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::getSlaveVertex( const Block &block ) const
{
    int32_t ctgid = block.getSlaveId();
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
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::PairedContigGraph( const std::list<Block> &blocks )
{
    this->initGraph( blocks );
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
void
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::initGraph( const std::list<Block> &blocks )
{
    this->initVertexLabels(blocks);

    // add an edge connected the master and slave contigs vertices of a block.
	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
        int32_t masterCtgId = b->getMasterId();
		int32_t slaveCtgId = b->getSlaveId();

        add_edge( this->_masterMap[masterCtgId], this->_slaveMap[slaveCtgId], *this );
    }
}


template <class VERTEX_PROP, class EDGE_WEIGHT>
void
PairedContigGraph<VERTEX_PROP,EDGE_WEIGHT>::initVertexLabels( const std::list<Block>& blocks )
{
    // collect master and slave contigs
    std::set< int32_t > masterCtgIdSet;
    std::set< int32_t > slaveCtgIdSet;

	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
		masterCtgIdSet.insert( b->getMasterId() );
		slaveCtgIdSet.insert( b->getSlaveId() );
    }

    this->_vertexToCtg.resize( masterCtgIdSet.size() + slaveCtgIdSet.size() );

    uint64_t i = 0;

    // insert into the vector the master contig IDs first.
	for( std::set< int32_t >::iterator id = masterCtgIdSet.begin(); id != masterCtgIdSet.end(); id++ )
    {
        this->_vertexToCtg.at(i) = *id;
        this->_masterMap[*id] = i;
        add_vertex(*this);
        i++;
    }

    this->_firstSlaveVertex = i; // index of the first slave contig vertex

    // insert into the vector the remaining contig IDs (of slave assembly).
    for( std::set< int32_t >::iterator id = slaveCtgIdSet.begin(); id != slaveCtgIdSet.end(); id++ )
    {
        this->_vertexToCtg.at(i) = *id;
        this->_slaveMap[*id] = i;
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
		os << "\t" << *v << "[label=\"" << this->_vertexToCtg[*v] << "\"];" << std::endl;

    EdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::edges(*this);

    for( EdgeIterator e = ebegin; e != eend; e++ )
		os << "\t" << boost::source(*e,*this) << "--" << boost::target(*e,*this) << "[color=black];" << std::endl;

    os << "}" << std::endl;

    return os;
}


#endif // PAIRED_GRAPH_CODE_