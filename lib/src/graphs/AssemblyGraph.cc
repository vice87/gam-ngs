/* 
 * File:   AssemblyGraph.cc
 * Author: vice
 * 
 * Created on 25 maggio 2011, 12.07
 */

#include "graphs/AssemblyGraph.hpp"

#include <vector>
#include <iostream>

#include <boost/graph/strong_components.hpp> 
#include <boost/graph/topological_sort.hpp> 
#include <boost/graph/graphviz.hpp>

#include "strand_fixer/RelativeStrand.hpp"
#include "OrderingFunctions.hpp"


AssemblyGraph::AssemblyGraph()
{
    std::vector< Block > blocks;
    this->initGraph( blocks );
}


AssemblyGraph::AssemblyGraph( const std::vector< Block > &blocks )
{
    this->initGraph( blocks );
}


const AssemblyGraph& 
AssemblyGraph::operator =(const AssemblyGraph& orig)
{
    *((Graph *)this) = *((Graph *)&orig);
    this->_blockVector = orig._blockVector;
    
    return *this;
}


const std::vector<Block>& 
AssemblyGraph::getBlocksVector() const
{
    return this->_blockVector;
}

const Block&
AssemblyGraph::getBlock(const UIntType& pos) const
{
    return this->_blockVector.at(pos);
}


std::list<Block>
AssemblyGraph::removeCyclesFromSCC( std::list<Block> &sccBlocks )
{
    if( sccBlocks.size() <= 1 ) return sccBlocks;
    
    int i = 0;
    
    // remove block with the minimum number of reads
    ReadNumBlocksOrderer rno;
    sccBlocks.erase( std::min_element( sccBlocks.begin(), sccBlocks.end(), rno) );
    
    // copy block list into a vector
    std::vector<Block> blockVect( sccBlocks.size() );
    for( std::list<Block>::iterator block = sccBlocks.begin(); block != sccBlocks.end(); block++ )
    {
        blockVect[i] = *block;
        i++;
    }    
    
    // build graph of assemblies
    AssemblyGraph graph( blockVect );
    
    typedef std::vector< size_t > container;
    container c;
    
    try
    {
        boost::topological_sort( graph, std::back_inserter(c) );
    }
    catch( boost::not_a_dag ) // if the graph still contains cycles
    {
        sccBlocks.clear();
        
        std::vector< UIntType > component( boost::num_vertices(graph) ),
                            discoverTime( boost::num_vertices(graph) ),
                            root( boost::num_vertices(graph) );
        std::vector< boost::default_color_type > color( boost::num_vertices(graph) );
        
        // compute strongly connected components
        UIntType sccNum = boost::strong_components( graph, &component[0],
                boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discoverTime[0]) );
        
        // group vertices by scc
        std::vector< std::list<Block> > components(sccNum);
        for( UIntType h=0; h != component.size(); ++h ) components[ component[h] ].push_back( blockVect[h] );
        
        // recursively remove cycles from SCCs
        for( UIntType h=0; h != sccNum; ++h )
        {            
            if( components[h].size() > 1 )
            {
                std::list<Block> tmpList = removeCyclesFromSCC( components[h] );
                sccBlocks.splice( sccBlocks.end(), tmpList );
            }
            else sccBlocks.splice( sccBlocks.end(), components[h] );
        }
    }
    
    return sccBlocks;    
}


void //std::vector< std::vector<Block> >
AssemblyGraph::removeCycles()
{
    std::vector< UIntType > component( boost::num_vertices(*this) ),
                            discoverTime( boost::num_vertices(*this) ),
                            root( boost::num_vertices(*this) );
    std::vector< boost::default_color_type > color( boost::num_vertices(*this) );
    
    // compute strongly connected components
    UIntType sccNum = boost::strong_components( *this, &component[0], 
            boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discoverTime[0]) );
    
    ////////////////////////////////////////////////////////// NEW CODE
    
    /*std::list< Block > newBlocks;
    
    // group vertices by scc
    std::vector< std::list<Block> > components(sccNum);
    for( UIntType h=0; h != component.size(); ++h )
        components.at( component.at(h) ).push_back( this->_blockVector.at(h) );
    
    for( UIntType h=0; h != sccNum; ++h )
    {
        if( components[h].size() > 1 )
        {
            std::list<Block> tmpList = removeCyclesFromSCC( components[h] );
            newBlocks.splice( newBlocks.end(), tmpList );
        }
        else newBlocks.splice( newBlocks.end(), components[h] );
    }
    
    UIntType i = 0;
    std::vector<Block> newBlockVect( newBlocks.size() );
    for( std::list<Block>::iterator block = newBlocks.begin(); block != newBlocks.end(); block++ )
    {
        newBlockVect[i] = *block;
        i++;
    }
    
    // re-initialize the graph with blocks that do not create cycles
    this->initGraph( newBlockVect );
    
    //std::cout << "Filtered blocks after new cycles removal = " << newBlockVect.size() << std::endl << std::flush;*/
    
    ////////////////////////////////////////////////////////// OLD CODE
    
    // group vertices by scc
    std::vector< std::list<UIntType> > components(sccNum);
    for( UIntType h=0; h != component.size(); ++h )
        components.at( component.at(h) ).push_back(h);
    
    // compute the number of strong components which include more than one vertex
    UIntType scc = 0;
    for( UIntType h=0; h != sccNum; ++h )
        if( components[h].size() > 1 ) scc++;
    
    std::vector< std::list<Block> > sComponents(scc);
    scc = 0;
    
    // fill strong components data structure
    UIntType newSize = this->_blockVector.size();
    for( UIntType h=0; h != sccNum; ++h )
    {
        if( components[h].size() > 1 )
        {
            newSize = newSize - components[h].size();
            
            std::list<UIntType>::iterator j;
            for( j = components[h].begin(); j != components[h].end(); j++ )
                sComponents[scc].push_back( this->_blockVector.at(*j) );
            
            scc++;
        }
    }
    
    // remove strong components from graph
    std::vector<Block> newBlockVector(newSize);
    UIntType newPos = 0;
    for( UIntType i = 0; i < this->_blockVector.size(); i++ )
    {
        if( components.at( component.at(i) ).size() == 1 )
        {
            newBlockVector.at(newPos) = this->_blockVector.at(i);
            newPos++;
        }
    }
    
    // re-initialize the graph with blocks that do not create cycles
    this->initGraph( newBlockVector );
    
    //std::cout << "Filtered blocks after cycles removal = " << newBlockVector.size() << std::endl << std::flush;
    
    //return sComponents;
}


void 
AssemblyGraph::initGraph( const std::vector<Block> &blocks )
{
    this->clear();
    this->_blockVector = blocks;
    
    StrandProbMap masterStrandMap, slaveStrandMap;
      
    boost::tie( masterStrandMap, slaveStrandMap ) =
            computeRelativeStrandMap( blocks );
    
    // index[i] -> index of i-th element in blocks
    // backIndex[i] -> index in the ordered vector of blocks[i]
    std::vector< UIntType > indexMaster, backIndexMaster;
    std::vector< UIntType > indexSlave, backIndexSlave;
    
    // compute master/slave blocks ordering and back ordering
    boost::tie( indexMaster, backIndexMaster ) = getOrderedMasterIndices( blocks );
    boost::tie( indexSlave, backIndexSlave ) = getOrderedSlaveIndices( blocks );
    
    // add a vertex for each block.
    this->addVertices();
    
    // for each block, connect his vertex to the successive blocks' vertices
    for( UIntType i=0; i < blocks.size(); i++ ) this->addMasterEdges( i, masterStrandMap, indexMaster, backIndexMaster );
    for( UIntType i=0; i < blocks.size(); i++ ) this->addSlaveEdges( i, slaveStrandMap, indexSlave, backIndexSlave );
}


UIntType 
AssemblyGraph::addVertices()
{
    std::vector<Block>::const_iterator block;
    for( block = this->_blockVector.begin(); block != this->_blockVector.end(); block++ )
    {
        Vertex v = boost::add_vertex(*this);
        // boost::put( boost::vertex_kind_t(), *this, v, both_vertex );
    }
    
    return this->_blockVector.size();
}


void
AssemblyGraph::addMasterEdges(
        const UIntType& vertex, 
        const StrandProbMap& strandMap, 
        const std::vector<UIntType>& index, 
        const std::vector<UIntType>& backIndex)
{
    IdType ctgId = (this->_blockVector[vertex]).getMasterFrame().getContigId();
    UIntType idx = backIndex.at(vertex);
    UIntType next = vertex, prev = vertex;
    
    if( idx+1 < index.size() ) next = index.at( idx+1 );
    if( idx > 0 ) prev = index.at( idx-1 );
    
    char strand = (strandMap.find(ctgId)->second).getStrand();
    
    // strand is '-' if it's more likely that the master contig is
    // reverse complemented, '+' otherwise
    switch( strand )
    {
        case '-': 
            std::swap(next,prev);
        case '+':
            if( next != vertex ) this->addMasterSingleEdge(vertex,next);
            if( prev != vertex ) this->addMasterSingleEdge(prev,vertex);
    }
}

void
AssemblyGraph::addSlaveEdges(
        const UIntType& vertex, 
        const StrandProbMap& strandMap, 
        const std::vector<UIntType>& index, 
        const std::vector<UIntType>& backIndex)
{
    IdType ctgId = (this->_blockVector[vertex]).getSlaveFrame().getContigId();
    UIntType idx = backIndex.at(vertex);
    UIntType next = vertex, prev = vertex;
    
    if( idx+1 < index.size() ) next = index.at( idx+1 );
    if( idx > 0 ) prev = index.at( idx-1 );
    
    char strand = (strandMap.find(ctgId)->second).getStrand();
    
    // strand is '-' if it's more likely that the slave contig is
    // reverse complemented, '+' otherwise
    switch( strand )
    {
        case '-': 
            std::swap(next,prev);
        case '+':
            if( next != vertex ) this->addSlaveSingleEdge(vertex,next);
            if( prev != vertex ) this->addSlaveSingleEdge(prev,vertex);
    }
}


bool
AssemblyGraph::addMasterSingleEdge(const UIntType& s, const UIntType& t)
{
    if( Block::shareMasterContig(this->_blockVector[s],this->_blockVector[t]) )
    {
        Edge e;
        bool exists;
        
        boost::tie(e,exists) = boost::edge(s,t,*this);
        if( !exists )
        {
            e = boost::add_edge(s,t,*this).first;
            boost::put( boost::edge_kind_t(), *this, e, MASTER_EDGE );
        }
        
        return true;
    }
    
    return false;
}

bool
AssemblyGraph::addSlaveSingleEdge(const UIntType& s, const UIntType& t)
{
    if( Block::shareSlaveContig(this->_blockVector[s],this->_blockVector[t]) )
    {
        Edge e;
        bool exists;
        
        boost::tie(e,exists) = boost::edge(s,t,*this);
        if( !exists ) 
        {
            e = boost::add_edge(s,t,*this).first;
            boost::put( boost::edge_kind_t(), *this, e, SLAVE_EDGE );
        }
        else
        {
            EdgeKindType edge_type = boost::get( boost::edge_kind_t(), *this, e );
            
            if( edge_type == MASTER_EDGE )
                boost::put( boost::edge_kind_t(), *this, e, BOTH_EDGE );
        }
        
        return true;
    }
    
    return false;
}


std::ostream&
AssemblyGraph::writeGraphviz(std::ostream& os)
{
    typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<AssemblyGraph>::edge_iterator EdgeIterator;
    
    os << "digraph AssemblyGraph {" << std::endl;
    os << "   rankdir=LR;" << std::endl;
    
    VertexIterator vbegin,vend;
    boost::tie(vbegin,vend) = boost::vertices(*this);
    
    for (VertexIterator v=vbegin; v!=vend; v++) 
    {
        IdType masterCtgId = this->_blockVector[*v].getMasterFrame().getContigId();
        UIntType masterFrameLen = this->_blockVector[*v].getMasterFrame().getLength();
        UIntType masterFrameBeg = this->_blockVector[*v].getMasterFrame().getBegin();
        
        IdType slaveCtgId = this->_blockVector[*v].getSlaveFrame().getContigId();
        UIntType slaveFrameLen = this->_blockVector[*v].getSlaveFrame().getLength();
        UIntType slaveFrameBeg = this->_blockVector[*v].getSlaveFrame().getBegin();
        
        os << "   " << *v << "[label=\"" 
                        << masterCtgId << ":" << masterFrameBeg << ":" << masterFrameLen << "\\n"
                        << slaveCtgId << ":" << slaveFrameBeg << ":" << slaveFrameLen 
                    << "\""
                    << ((boost::in_degree(*v,*this) > 1 || boost::out_degree(*v,*this) > 1) ? ", color = blue" : "") << "];"
                    << std::endl;
    }
    
    EdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::edges(*this);
    
    for( EdgeIterator e = ebegin; e != eend; e++ )
    {
        EdgeKindType kind = boost::get(boost::edge_kind_t(), *this, *e);
        
        switch(kind)
        {
            case MASTER_EDGE:
                os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=black];" << std::endl;
                break;
            case SLAVE_EDGE:
                os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=red];" << std::endl;
                break;
            case BOTH_EDGE:
                os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=green];" << std::endl;
                break;
        }
    }
    
    os << "}" << std::endl;
    
    return os;
}


void AssemblyGraph::agTopologicalSort( const AssemblyGraph &g, Vertex v, std::vector<char> &colors, std::list<Vertex> &tsList )
{
    colors[v] = 1;
    
    AdjacencyIterator begin, end;
    boost::tie(begin,end) = boost::adjacent_vertices(v,g);
    
    for( AdjacencyIterator u = begin; u != end; u++ )
    {
        if( colors[*u] == 0 ) AssemblyGraph::agTopologicalSort( g, *u, colors, tsList );
        if( colors[*u] == 1 ) throw boost::not_a_dag();
    }
    
    colors[v] = 2;
    tsList.push_back(v);
}


void AssemblyGraph::agTopologicalSort( const AssemblyGraph &g, std::list<Vertex> &tsList )
{    
    std::list<Vertex> roots;
    std::vector<char> colors( boost::num_vertices(g), 0 );
    
    VertexIterator vbegin,vend;
    boost::tie(vbegin,vend) = boost::vertices(g);
    
    for( VertexIterator v = vbegin; v != vend; v++ ) 
    {
        if( boost::in_degree(*v,g) == 0 ) roots.push_back(*v);
    }
    
    for( std::list<Vertex>::const_iterator v = roots.begin(); v != roots.end(); v++ )
    {
        AssemblyGraph::agTopologicalSort(g, *v, colors, tsList);
    }
}