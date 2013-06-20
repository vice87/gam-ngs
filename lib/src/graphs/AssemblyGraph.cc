/*
 * File:   AssemblyGraph.cc
 * Author: vice
 *
 * Created on 25 maggio 2011, 12.07
 */

#include "graphs/AssemblyGraph.hpp"

AssemblyGraph::AssemblyGraph( uint64_t id ) : _agId(id)
{	
    std::list< Block > blocks;
    this->initGraph( blocks );
}


AssemblyGraph::AssemblyGraph( const std::list< Block > &blocks, uint64_t id ) : _agId(id)
{
    this->initGraph( blocks );
}


const AssemblyGraph&
AssemblyGraph::operator =(const AssemblyGraph& orig)
{
    *((Graph *)this) = *((Graph *)&orig);
    this->_blockVector = orig._blockVector;
	this->_agId = orig._agId;

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


void
AssemblyGraph::removeForks()
{
    // remove nodes with in-degree greater than 1
    std::set< std::pair<IdType,IdType> > badBlocks;
	std::vector< Vertex > fork_blocks;
    bool has_forks = false;

    VertexIterator vbegin,vend;
    boost::tie(vbegin,vend) = boost::vertices(*this);

	for( VertexIterator v = vbegin; v != vend; v++ )
    {
        if( boost::in_degree(*v,*this) > 1 || boost::out_degree(*v,*this) > 1 )
        {
            //IdType mID = this->_blockVector[*v].getMasterFrame().getContigId();
            //IdType sID = this->_blockVector[*v].getSlaveFrame().getContigId();
            //badBlocks.insert( std::make_pair(mID,sID) );
			fork_blocks.push_back(*v);
            has_forks = true;
        }
    }

    double min_cov;
	Vertex del_vtx;

    for( size_t i=0; i < fork_blocks.size(); i++ )
	{
		Frame &mf = (this->_blockVector[ fork_blocks[i] ]).getMasterFrame();
		Frame &sf = (this->_blockVector[ fork_blocks[i] ]).getSlaveFrame();

		if( i==0 )
		{
			min_cov = std::min( double(mf.getBlockReadsLen())/double(mf.getReadsLen()), double(sf.getBlockReadsLen())/double(sf.getReadsLen()) );
			del_vtx = fork_blocks[0];
		}
		else
		{
			double cur_cov = std::min( double(mf.getBlockReadsLen())/double(mf.getReadsLen()), double(sf.getBlockReadsLen())/double(sf.getReadsLen()) );
			if( cur_cov < min_cov )
			{
				min_cov = cur_cov;
				del_vtx = fork_blocks[i];
			}
		}
	}

	std::list<Block> newBlocks;

	for( VertexIterator v = vbegin; v != vend; v++ )
	{
		if( *v != del_vtx ) newBlocks.push_back( this->_blockVector[*v] );
	}

	this->initGraph( newBlocks );
	if(has_forks) this->removeForks();

	return;

    /*for( VertexIterator v = vbegin; v != vend; v++ )
    {
        IdType mID = this->_blockVector[*v].getMasterFrame().getContigId();
        IdType sID = this->_blockVector[*v].getSlaveFrame().getContigId();

        if( badBlocks.find( std::make_pair(mID,sID) ) == badBlocks.end() )
            newBlocks.push_back( this->_blockVector[*v] );
    }

    this->initGraph( newBlocks );

    if(has_forks) this->removeForks();*/
}


void
AssemblyGraph::initGraph( const std::list<Block> &blocks )
{
    this->clear();

	// fill local blocks vector
	this->_blockVector.reserve( blocks.size() );
	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); ++b ) this->_blockVector.push_back(*b);

    StrandProbMap masterStrandMap, slaveStrandMap;

    boost::tie( masterStrandMap, slaveStrandMap ) = computeRelativeStrandMap( blocks );

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
    int32_t ctgId = (this->_blockVector[vertex]).getMasterId();
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
    int32_t ctgId = (this->_blockVector[vertex]).getSlaveId();
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
			EdgeProperty edge_prop = { MASTER_EDGE, 0.0 };
			e = boost::add_edge(s,t,*this).first;
			boost::put( boost::edge_kind_t(), *this, e, edge_prop );
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
			EdgeProperty edge_prop = { SLAVE_EDGE, 0.0 };
            e = boost::add_edge(s,t,*this).first;
			boost::put( boost::edge_kind_t(), *this, e, edge_prop );
        }
        else
        {
            //EdgeKindType edge_type = boost::get( boost::edge_kind_t(), *this, e );
            EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), *this, e );
			EdgeKindType edge_type = edge_prop.kind;

            if( edge_type == MASTER_EDGE )
			{
				edge_prop.kind = BOTH_EDGE;
                boost::put( boost::edge_kind_t(), *this, e, edge_prop );
			}
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
        IdType masterAId = this->_blockVector[*v].getMasterFrame().getAssemblyId();
        IdType masterCtgId = this->_blockVector[*v].getMasterFrame().getContigId();
        UIntType masterFrameLen = this->_blockVector[*v].getMasterFrame().getLength();
        UIntType masterFrameBeg = this->_blockVector[*v].getMasterFrame().getBegin();
        RealType mbc = ((RealType) this->_blockVector[*v].getMasterFrame().getBlockReadsLen()) / ((RealType)masterFrameLen);
        RealType mc = ((RealType) this->_blockVector[*v].getMasterFrame().getReadsLen()) / ((RealType)masterFrameLen);

        IdType slaveAId = this->_blockVector[*v].getSlaveFrame().getAssemblyId();
        IdType slaveCtgId = this->_blockVector[*v].getSlaveFrame().getContigId();
        UIntType slaveFrameLen = this->_blockVector[*v].getSlaveFrame().getLength();
        UIntType slaveFrameBeg = this->_blockVector[*v].getSlaveFrame().getBegin();
        RealType sbc = ((RealType) this->_blockVector[*v].getSlaveFrame().getBlockReadsLen()) / ((RealType)slaveFrameLen);
        RealType sc = ((RealType) this->_blockVector[*v].getSlaveFrame().getReadsLen()) / ((RealType)slaveFrameLen);

        os << "   " << *v << "[label=\""
                    << "<" << masterAId << "," << masterCtgId << "> :" << masterFrameBeg << ":" << masterFrameLen << " (" << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mbc << "/" << mc << ")" << "\\n"
                    << "<" << slaveAId << "," << slaveCtgId << "> :" << slaveFrameBeg << ":" << slaveFrameLen << " (" << std::setiosflags(std::ios::fixed) << std::setprecision(2) << sbc << "/" << sc << ")"
                    << "\""
					<< ((boost::in_degree(*v,*this) > 1 || boost::out_degree(*v,*this) > 1) ? ", color=deepskyblue, style=filled" : "") << "];"
                    << std::endl;
    }

    EdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::edges(*this);

    for( EdgeIterator e = ebegin; e != eend; e++ )
    {
        //EdgeKindType kind = boost::get(boost::edge_kind_t(), *this, *e);
        EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), *this, *e);
		EdgeKindType kind = edge_prop.kind;
		double weight = edge_prop.weight;

        switch(kind)
        {
            case MASTER_EDGE:
				os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=black, label=\"" << weight << "\"];" << std::endl;
                break;
            case SLAVE_EDGE:
				os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=red, label=\"" << weight << "\"];" << std::endl;
                break;
            case BOTH_EDGE:
				os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=green, label=\"" << weight << "\"];" << std::endl;
                break;
        }
    }

    os << "}" << std::endl;

    return os;
}


void AssemblyGraph::reverseEdges()
{
    std::vector< std::pair<Vertex,Vertex> > new_edges;
    std::vector< EdgeProperty > edge_properties; //std::vector< EdgeKindType > types;

    EdgeIterator begin,end;
    boost::tie(begin,end) = boost::edges(*this);
    for( EdgeIterator e = begin; e != end; e++ )
    {
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), *this, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), *this, *e);
        new_edges.push_back( std::make_pair( boost::target(*e,*this), boost::source(*e,*this) ) );
        edge_properties.push_back(edge_prop); //types.push_back( kind );
    }

    size_t vertices = boost::num_vertices(*this);

    this->clear();
    for( size_t i=0; i < vertices; i++ ) boost::add_vertex(*this);

    for( size_t i=0; i < new_edges.size(); i++ )
    {
        Edge e = boost::add_edge( new_edges[i].first, new_edges[i].second, *this ).first;
		boost::put( boost::edge_kind_t(), *this, e, edge_properties[i] );
    }
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


bool AssemblyGraph::hasForks()
{
	typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
	typedef boost::graph_traits<AssemblyGraph>::edge_iterator EdgeIterator;

	VertexIterator vbegin,vend;
	boost::tie(vbegin,vend) = boost::vertices(*this);

	for (VertexIterator v=vbegin; v!=vend; v++)
	{
		int in_deg = boost::in_degree(*v,*this);
		int out_deg = boost::out_degree(*v,*this);
		
		if( in_deg > 1 || out_deg > 1 ) return true;
	}
	
	return false;
}


bool AssemblyGraph::hasBubbles()
{
	std::list<Vertex> roots;
	std::vector<char> colors( boost::num_vertices(*this), 0 );
	bool found = false;

	VertexIterator vbegin,vend;
	boost::tie(vbegin,vend) = boost::vertices(*this);

	// find roots
	for( VertexIterator v = vbegin; v != vend; v++ ) if( boost::in_degree(*v,*this) == 0 )
		roots.push_back(*v);

	for( std::list<Vertex>::const_iterator v = roots.begin(); v != roots.end(); v++ )
	{
		for(size_t i = 0; i < colors.size(); i++) colors[i] = 0;
		this->bubbleDFS(*v, colors, found);
	}

	return found;
}

void AssemblyGraph::bubbleDFS( Vertex v, std::vector<char> &colors, bool &found )
{
	colors[v] = 1;

	AdjacencyIterator begin, end;
	boost::tie(begin,end) = boost::adjacent_vertices(v,*this);

	for( AdjacencyIterator u = begin; u != end; u++ )
	{
		if( colors[*u] == 0 ) this->bubbleDFS( *u, colors, found );
		else if( colors[*u] == 2 ) found = true;
		else if( colors[*u] == 1 ) throw boost::not_a_dag();
	}

	colors[v] = 2;
}