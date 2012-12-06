/*
 * File:   AssemblyGraph.cc
 * Author: vice
 *
 * Created on 25 maggio 2011, 12.07
 */

#include "graphs/CompactAssemblyGraph.hpp"

#include <vector>
#include <iostream>
#include <iomanip>

#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>

#include "strand_fixer/RelativeStrand.hpp"
#include "OrderingFunctions.hpp"


CompactAssemblyGraph::CompactAssemblyGraph( const AssemblyGraph &ag )
{
    this->initGraph(ag);
}


const CompactAssemblyGraph&
CompactAssemblyGraph::operator =(const CompactAssemblyGraph& orig)
{
    *((Graph *)this) = *((Graph *)&orig);
    this->_blockVector = orig._blockVector;

    return *this;
}


const std::vector< std::list<Block> >&
CompactAssemblyGraph::getBlocksVector() const
{
    return this->_blockVector;
}


const std::list<Block>&
CompactAssemblyGraph::getBlocks( const Vertex &pos ) const
{
    return this->_blockVector.at(pos);
}


void
CompactAssemblyGraph::initGraphDFS( const AssemblyGraph &ag, Vertex v, std::vector<char> &colors, std::vector<Vertex> &ag2cg, Vertex u )
{
    if( colors[v] == 1 ) // if node already visited, add edge and return
    {
        Edge e;
        bool exists;
        boost::tie(e,exists) = boost::edge(u,v,ag);
		EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), ag, e ); //EdgeKindType edge_type = boost::get( boost::edge_kind_t(), ag, e );

        e = boost::add_edge( ag2cg[u], ag2cg[v], *this ).first;
		boost::put( boost::edge_kind_t(), *this, e, edge_prop );

        return;
    }

	colors[v] = 1; // mark current node as visited

	Edge e;
	bool exists;
	boost::tie(e,exists) = boost::edge(u,v,ag);

	EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), ag, e ); //EdgeKindType edge_type = boost::get( boost::edge_kind_t(), ag, e );

	if( edge_prop.kind == BOTH_EDGE ) // if previous vertex is connected to the current one with a BOTH_EDGE (master+slave edge)
	{
		_blockVector.at( ag2cg[u] ).push_back( ag.getBlock(v) );
		ag2cg[v] = ag2cg[u];
	}
	else // else add a new vertex to the compact graph
	{
		boost::add_vertex(*this);
		ag2cg[v] = _num_vertices;
		_blockVector.push_back( std::list<Block>(1,ag.getBlock(v)) );

		e = boost::add_edge( ag2cg[u], ag2cg[v], *this ).first;
		boost::put( boost::edge_kind_t(), *this, e, edge_prop );

		_num_vertices++;
	}

	AdjacencyIterator begin, end;
	boost::tie(begin,end) = boost::adjacent_vertices(v,ag);
	for( AdjacencyIterator z = begin; z != end; z++ ) this->initGraphDFS(ag, *z, colors, ag2cg, v);
}


void
CompactAssemblyGraph::initGraph( const AssemblyGraph &ag )
{
    this->clear();
	this->_num_vertices = 0;

	size_t ag_vertices = boost::num_vertices(ag);

	std::vector<char> colors( ag_vertices, 0 );
	std::vector<Vertex> ag2cg( ag_vertices, 0 ); // associate each vertex of ag to a vertex of this CompactAssemblyGraph

	VertexIterator vbegin,vend;
	boost::tie(vbegin,vend) = boost::vertices(ag);

	for( VertexIterator v = vbegin; v != vend; v++ )
	{
		if( boost::in_degree(*v,ag) == 0 )
        {
            boost::add_vertex(*this);
            colors[*v] = 1;
            ag2cg[*v] = _num_vertices;
            _blockVector.push_back( std::list<Block>(1,ag.getBlock(*v)) );
            _num_vertices++;

            AdjacencyIterator begin, end;
            boost::tie(begin,end) = boost::adjacent_vertices(*v,ag);
            for( AdjacencyIterator z = begin; z != end; z++ ) this->initGraphDFS(ag, *z, colors, ag2cg, *v);
        }
	}
}


void
CompactAssemblyGraph::computeEdgeWeights( MultiBamReader &masterBamReader, MultiBamReader &masterMpBamReader, 
										  MultiBamReader &slaveBamReader, MultiBamReader &slaveMpBamReader )
{
	EdgeIterator ebegin,eend;
	boost::tie(ebegin,eend) = boost::edges(*this);
	
	for( EdgeIterator e = ebegin; e != eend; e++ )
	{
		//EdgeKindType kind = boost::get(boost::edge_kind_t(), *this, *e);
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), *this, *e);
		EdgeKindType kind = edge_prop.kind;
		
		std::list<Block>& b1 = _blockVector.at( boost::source(*e,*this) );
		std::list<Block>& b2 = _blockVector.at( boost::target(*e,*this) );
		
		std::pair<double,int32_t> edge_lab;
		
		switch(kind)
		{
			case MASTER_EDGE:
				edge_lab = this->getRegionScore( masterBamReader, masterMpBamReader, MASTER_EDGE, b1, b2 );
				edge_prop.weight = edge_lab.first;
				edge_prop.rnum = edge_lab.second;
				break;
				
			case SLAVE_EDGE:
				edge_lab = this->getRegionScore( slaveBamReader, slaveMpBamReader, SLAVE_EDGE, b1, b2 );
				edge_prop.weight = edge_lab.first;
				edge_prop.rnum = edge_lab.second;
				break;
				
			default:
				edge_prop.weight = 0.0;
				edge_prop.rnum = 0;
				break;
		}
		
		// put edge weight
		boost::put( boost::edge_kind_t(), *this, *e, edge_prop );
	}
}


std::pair<double,int32_t>
CompactAssemblyGraph::getRegionScore( MultiBamReader &peBamReader, MultiBamReader &mpBamReader, EdgeKindType kind, 
									  std::list<Block>& b1, std::list<Block>& b2 )
{
	//std::cerr << "PE" << std::endl;
	//std::vector<double> peScore = getLibRegionScore( peBamReader, kind, b1, b2 );
	//std::cerr << "MP" << std::endl;
	//std::vector<double> mpScore;
	std::pair< std::vector<double>, std::vector<int32_t> > mpStats, peStats;
	if(mpBamReader.size() > 0) mpStats = getLibRegionScore2( mpBamReader, kind, b1, b2 );
	if(peBamReader.size() > 0) peStats = getLibRegionScore2( peBamReader, kind, b1, b2 );
	
	std::pair<double,int32_t> mp_edge_weight = std::make_pair(-10.0,0);
	std::pair<double,int32_t> pe_edge_weight = std::make_pair(-10.0,0);
	
	std::vector<double>  &mpWeights = mpStats.first;
	std::vector<int32_t> &mpNReads = mpStats.second;
	
	std::vector<double>  &peWeights = peStats.first;
	std::vector<int32_t> &peNReads = peStats.second;
	
	for( size_t i=0; i < mpWeights.size(); i++ )
	{
		if( mpNReads[i] > mp_edge_weight.second ) mp_edge_weight = std::make_pair(mpWeights[i],mpNReads[i]);
	}
	
	for( size_t i=0; i < peWeights.size(); i++ )
	{
		if( peNReads[i] > pe_edge_weight.second ) pe_edge_weight = std::make_pair(peWeights[i],peNReads[i]);
	}
	
	if( mp_edge_weight.second > 0 && pe_edge_weight.second == 0 ) return mp_edge_weight;
	if( pe_edge_weight.second > 0 && mp_edge_weight.second == 0 ) return pe_edge_weight;
	
	return (mp_edge_weight.first < pe_edge_weight.first) ? mp_edge_weight : pe_edge_weight;
	
	//return (mpStats.first.size() > 0) ? std::make_pair(mpStats.first.at(0),mpStats.second.at(0)) : std::make_pair(-10.0,0);
	
	//double maxScore = 0.0;
	//for( size_t i=0; i < peScore.size(); i++ ) if( peScore[i] > maxScore ) maxScore = peScore[i];
	//for( size_t i=0; i < mpScore.size(); i++ ) if( mpScore[i] > maxScore ) maxScore = mpScore[i];
	
	//return maxScore;
}


std::pair< std::vector<double>, std::vector<int32_t> >
CompactAssemblyGraph::getLibRegionScore2( MultiBamReader &bamReader, EdgeKindType kind, 
										  std::list<Block>& b1, std::list<Block>& b2 )
{
	int32_t id, seq_len, start, end, gap, region, s1, s2, t;
	uint64_t good_reads, exp_reads, num_reads;
	
	std::vector<double> score( bamReader.size(), -4 );
	std::vector<int32_t> rnum( bamReader.size(), 0 );
	const RefVector& ref = bamReader.GetReferenceData();
	
	// this shouldn't happen
	if( kind != MASTER_EDGE && kind != SLAVE_EDGE ) return std::make_pair(score,rnum);
	if( b1.size() == 0 || b2.size() == 0 ) return std::make_pair(score,rnum);
	
	// COMPUTE STATISTICS FOR EACH LIBRARY
	for( int lib=0; lib < bamReader.size(); lib++ )
	{
		int32_t isizeLibMean = bamReader.getISizeMean(lib);
		int32_t isizeLibStd = bamReader.getISizeStd(lib);
		
		int32_t minInsert = isizeLibMean - 3*isizeLibStd;
		int32_t maxInsert = isizeLibMean + 3*isizeLibStd;
		
		if(minInsert < 0) minInsert = 0;
		
		Frame& f1 = (kind == MASTER_EDGE) ? b1.front().getMasterFrame() : b1.front().getSlaveFrame();
		Frame& f2 = (kind == MASTER_EDGE) ? b2.front().getMasterFrame() : b2.front().getSlaveFrame();
		Frame& l1 = (kind == MASTER_EDGE) ? b1.back().getMasterFrame() : b1.back().getSlaveFrame();
		Frame& l2 = (kind == MASTER_EDGE) ? b2.back().getMasterFrame() : b2.back().getSlaveFrame();
		
		id = f1.getContigId();
		seq_len = ref[id].RefLength;
		
		int32_t r1_beg = std::min( f1.getBegin(), l1.getBegin() );
		int32_t r1_end = std::max( f1.getEnd(), l1.getEnd() );
		int32_t r2_beg = std::min( f2.getBegin(), l2.getBegin() );
		int32_t r2_end = std::max( f2.getEnd(), l2.getEnd() );
		
		// skip included frames
		if( (r1_beg <= r2_beg && r1_end >= r2_end) || 
			(r2_beg <= r1_beg && r2_end >= r1_end) )
		{
			score[lib] = -1;
			continue;
		}
		
		gap = (r1_beg <= r2_beg) ? (r2_beg - r1_end + 1) : (r1_beg - r2_end + 1);
		
		t = (r1_beg <= r2_beg) ? (gap >= 0 ? r2_beg : r1_end) : (gap >= 0 ? r1_beg : r2_end);
		s1 = std::max( t - maxInsert, 0 );
		s2 = (r1_beg <= r2_beg) ? (gap >= 0 ? r1_end : r2_beg) : (gap >= 0 ? r2_end : r1_beg);
		
		if( seq_len - s1 < maxInsert )
		{
			score[lib] = -2;
			continue;
		}
		
		if( gap >= maxInsert )
		{
			score[lib] = -3;
			continue;
		}
		
		//start = (f1_beg <= f2_beg) ? std::max( f1_end - maxInsert, 0 ) : std::max( f2_end - maxInsert, 0 ); //std::min( f1.getBegin(), f2.getBegin() );
		//end = std::max( f1_end, f2_end );
		//region = end - start + 1;
		
		//s1 = start; //(f1.getBegin() <= f2.getBegin()) ? std::max( f1.getEnd() - maxInsert, 0 ) : std::max( f2.getEnd() - maxInsert, 0 );
		//s2 = (f1_beg <= f2_beg) ? f1_end : f2_end;
		//t = (f1_beg <= f2_beg) ? f2_beg : f1_beg;
		
		// retrieve BAM readers for current library
		bamReader.lockBamReader(lib);
		
		BamReader *reader = bamReader.getBamReader(lib);
		reader->SetRegion( id, s1, id, s2+1 );
		
		int32_t nh, xt;
		good_reads = 0;
		exp_reads = 0;
		num_reads = 0;
		
		BamAlignment align;
		while( reader->GetNextAlignmentCore(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;
			//if( !align.IsMateMapped() || align.RefID != align.MateRefID || align.MatePosition < t ) continue;
			
			align.BuildCharData(); // fill string fields
			
			// if not defined, I assume read's multiplicity is 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard field
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';	// bwa field
			if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1
			
			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;
			int32_t startMate = align.MatePosition;
			int32_t endMate = startMate + readLength - 1;
			
			// don't count reads not completely included in the region
			if( startRead < s1 || startRead > s2 ) continue; //|| endRead > s2 ) continue;
			
			if( !align.IsReverseStrand() )
			{
				int32_t minInsertPos = startRead + minInsert;
				int32_t maxInsertPos = startRead + maxInsert;
				int32_t readOverlap = endRead > s2 ? s2-startRead+1 : readLength;
				
				/*if( !align.IsMateMapped() ){ reads++; continue; }
				if( align.RefID != align.MateRefID ){ if(endPos2 < seq_len) reads++; continue; }
				if( align.IsMateReverseStrand() && startMate >= t ){ good_reads++; reads++; }*/
				
				// unmapped mate
				if( !align.IsMateMapped() ){ exp_reads += readOverlap; num_reads++; continue; }
				// mate is mapped in a different sequence while it should not be.
				if( align.RefID != align.MateRefID ){ if( maxInsertPos < seq_len ) exp_reads += readOverlap; num_reads++; continue; }
				// mate mapped in the same sequence, crossing the gap, with wrong orientation
				if( !align.IsMateReverseStrand() && endMate >= t ){ exp_reads += readOverlap; num_reads++;}
				// mate mapped in the same sequence, crossing the gap, with correct orientation
				if( align.IsMateReverseStrand() && endMate >= t ){ good_reads += readOverlap; exp_reads += readOverlap; num_reads++; }
			}
		} // end while
		
		bamReader.unlockBamReader(lib);
		
		if( num_reads < 10 || exp_reads == 0 )
		{ 
			score[lib] = -5; 
			rnum[lib] = 0;
		}
		else
		{
			score[lib] = good_reads / ((double)exp_reads);
			rnum[lib] = num_reads;
		}
	}
	
	return std::make_pair(score,rnum);
}


bool CompactAssemblyGraph::hasBubbles()
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

bool CompactAssemblyGraph::bubbleDFS( Vertex v, std::vector<char> &colors, bool &found )
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



std::ostream&
CompactAssemblyGraph::writeGraphviz(std::ostream& os)
{
    os << "digraph CompactAssemblyGraph {" << std::endl;
    os << "   rankdir=LR;" << std::endl;

    VertexIterator vbegin,vend;
    boost::tie(vbegin,vend) = boost::vertices(*this);

    for (VertexIterator v=vbegin; v!=vend; v++)
    {
        const std::list<Block>& listRef = this->_blockVector[*v];
        const Block& firstBlock = listRef.front();
        const Block& lastBlock = listRef.back();

        int32_t masterId = firstBlock.getMasterId();
        int32_t slaveId = firstBlock.getSlaveId();

        int32_t masterBegin = std::min( firstBlock.getMasterFrame().getBegin(), lastBlock.getMasterFrame().getBegin() );
        int32_t masterEnd = std::max( firstBlock.getMasterFrame().getEnd(), lastBlock.getMasterFrame().getEnd() );

        int32_t slaveBegin = std::min( firstBlock.getSlaveFrame().getBegin(), lastBlock.getSlaveFrame().getBegin() );
        int32_t slaveEnd = std::max( firstBlock.getSlaveFrame().getEnd(), lastBlock.getSlaveFrame().getEnd() );

        size_t blocks = listRef.size();

        os << "   " << *v << "[label=\""
                    << masterId << ": " << masterBegin << "-" << masterEnd << "\\n"
                    << slaveId << ": " << slaveBegin << "-" << slaveEnd << "\\n"
                    << "blocks = " << blocks
                    << "\""
                    << ((boost::in_degree(*v,*this) > 1 || boost::out_degree(*v,*this) > 1) ? ", color=deepskyblue, style=filled " : "") << "];"
                    << std::endl;
    }

    EdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::edges(*this);

    for( EdgeIterator e = ebegin; e != eend; e++ )
    {
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), *this, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), *this, *e);

		EdgeKindType kind = edge_prop.kind;
		double weight = edge_prop.weight;
		int32_t rnum = edge_prop.rnum;

        switch(kind)
        {
            case MASTER_EDGE:
				os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=black, label=\"" << weight << "/" << rnum << "\"];" << std::endl;
                break;
            case SLAVE_EDGE:
				os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=red, label=\"" << weight << "/" << rnum << "\"];" << std::endl;
                break;
            case BOTH_EDGE:
				os << "   " << boost::source(*e,*this) << "->" << boost::target(*e,*this) << "[color=green, label=\"" << weight << "\"];" << std::endl;
                break;
        }
    }

    os << "}" << std::endl;

    return os;
}