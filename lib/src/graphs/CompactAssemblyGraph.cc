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

#include "graphs/CompactAssemblyGraph.hpp"

#include <stack>

#include "OptionsMerge.hpp"
using namespace options;
extern OptionsMerge g_options;


CompactAssemblyGraph::CompactAssemblyGraph( const AssemblyGraph &ag )
{
	this->_cgId = ag.getId();
    this->initGraph2(ag);
}


const CompactAssemblyGraph&
CompactAssemblyGraph::operator =(const CompactAssemblyGraph& orig)
{
    *((CompactAssemblyGraph *)this) = *((CompactAssemblyGraph *)&orig);
    this->_blockVector = orig._blockVector;
	this->_cgId = orig._cgId;
    this->_num_vertices = orig._num_vertices;

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
CompactAssemblyGraph::initGraphDFS_NR(
	const AssemblyGraph &ag,
	const AssemblyGraph::Vertex &root,
	boost::dynamic_bitset<> *visited,
	std::vector<Vertex> *ag2cg )
{
	AssemblyGraph::Edge e_ag;
	AssemblyGraph::AdjacencyIterator begin, end;
	Edge e; bool exists; Vertex new_v;

	std::stack<AssemblyGraph::Vertex> *cur_stack = new std::stack<AssemblyGraph::Vertex>();
	std::stack<AssemblyGraph::Vertex> *pre_stack = new std::stack<AssemblyGraph::Vertex>();

	new_v = boost::add_vertex(*this);
	std::list<Block> new_list = std::list<Block>(1,ag.getBlock(root));
	_blockVector.push_back( new_list );

	visited->set(root);
	ag2cg->at(root) = new_v;
	this->_num_vertices++;

	boost::tie(begin,end) = boost::adjacent_vertices(root,ag);
	for( AssemblyGraph::AdjacencyIterator z = begin; z != end; z++ )
	{
		cur_stack->push(*z);
		pre_stack->push(root);
	}

	while( not cur_stack->empty() )
	{
		AssemblyGraph::Vertex curr = cur_stack->top(); cur_stack->pop();
		AssemblyGraph::Vertex prev = pre_stack->top(); pre_stack->pop();

		if( visited->test(curr) ) // if node already visited, add edge and continue
		{
			boost::tie(e_ag,exists) = boost::edge(prev,curr,ag);
			EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), ag, e_ag );

			boost::tie(e,exists) = boost::add_edge( ag2cg->at(prev), ag2cg->at(curr), *this );
			boost::put( boost::edge_kind_t(), *this, e, edge_prop );

			continue;
		}

		// if curr has not been visited yet
		visited->set(curr); // mark current node as visited

		boost::tie(e_ag,exists) = boost::edge( prev, curr, ag );
		EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), ag, e_ag );

		if( edge_prop.kind == BOTH_EDGE ) // if previous vertex is connected to the current one with a BOTH_EDGE (master+slave edge)
		{
			_blockVector.at( ag2cg->at(prev) ).push_back( ag.getBlock(curr) );
			ag2cg->at(curr) = ag2cg->at(prev);
		}
		else // else add a new vertex to the compact graph
		{
			Vertex new_v = boost::add_vertex(*this);
			ag2cg->at(curr) = new_v;
			_blockVector.push_back( std::list<Block>(1,ag.getBlock(curr)) );

			this->_num_vertices++;

			//std::cerr << "debug: vertices=" << this->_num_vertices << "\n";
			//std::cerr << "debug: curr=" << curr << " ag2cg(curr)=" << ag2cg->at(curr) << "\n";
			//std::cerr << "debug: prev=" << prev << " ag2cg(prev)=" << ag2cg->at(prev) << std::endl;

			boost::tie(e,exists) = boost::add_edge( ag2cg->at(prev), ag2cg->at(curr), *this );
			boost::put( boost::edge_kind_t(), *this, e, edge_prop );
		}

		boost::tie(begin,end) = boost::adjacent_vertices(curr,ag);
		for( AssemblyGraph::AdjacencyIterator z = begin; z != end; z++ )
		{
			cur_stack->push(*z);
			pre_stack->push(curr);
		}
	}

	delete cur_stack;
	delete pre_stack;
}

void
CompactAssemblyGraph::initGraph2( const AssemblyGraph &ag )
{
	this->clear();
	this->_num_vertices = 0;

	size_t ag_vertices = boost::num_vertices(ag);

	boost::dynamic_bitset<> *visited = new boost::dynamic_bitset<>(ag_vertices);
	std::vector<Vertex> *ag2cg = new std::vector<Vertex>(ag_vertices,0);

	AssemblyGraph::VertexIterator vbegin,vend;
	boost::tie(vbegin,vend) = boost::vertices(ag);

	for( AssemblyGraph::VertexIterator r = vbegin; r != vend; r++ )
	{
        if( boost::in_degree(*r,ag) == 0 && !visited->test(*r) ) // for each unvisited root
        {
            this->initGraphDFS_NR( ag, *r, visited, ag2cg );
        }
	}

	delete visited;
	delete ag2cg;
}


void
CompactAssemblyGraph::initGraphDFS(
	const AssemblyGraph &ag,
	const AssemblyGraph::Vertex &v,
	const AssemblyGraph::Vertex &u,
	boost::dynamic_bitset<> *colors,
	std::vector<Vertex> *ag2cg )
{
	Edge e; bool exists;

    if( colors->test(v) ) // if node already visited, add edge and return
    {
        boost::tie(e,exists) = boost::edge(u,v,ag);
		EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), ag, e ); //EdgeKindType edge_type = boost::get( boost::edge_kind_t(), ag, e );

        e = boost::add_edge( ag2cg->at(u), ag2cg->at(v), *this ).first;
		boost::put( boost::edge_kind_t(), *this, e, edge_prop );

        return;
    }

	colors->set(v); // mark current node as visited

	boost::tie(e,exists) = boost::edge(u,v,ag);
	EdgeProperty edge_prop = boost::get( boost::edge_kind_t(), ag, e ); //EdgeKindType edge_type = boost::get( boost::edge_kind_t(), ag, e );

	if( edge_prop.kind == BOTH_EDGE ) // if previous vertex is connected to the current one with a BOTH_EDGE (master+slave edge)
	{
		_blockVector.at( ag2cg->at(u) ).push_back( ag.getBlock(v) );
		ag2cg->at(v) = ag2cg->at(u);
	}
	else // else add a new vertex to the compact graph
	{
		Vertex new_v = boost::add_vertex(*this);
		ag2cg->at(v) = new_v;
		_blockVector.push_back( std::list<Block>(1,ag.getBlock(v)) );

		_num_vertices++;

		std::cerr << "warning: u vertex-desc=" << ag2cg->at(u) << " num-vertices=" << _num_vertices << std::endl;
		std::cerr << "warning: v vertex-desc=" << ag2cg->at(v) << " num-vertices=" << _num_vertices << std::endl;

		boost::tie(e,exists) = boost::add_edge( ag2cg->at(u), ag2cg->at(v), *this );
		boost::put( boost::edge_kind_t(), *this, e, edge_prop );
	}

	AssemblyGraph::AdjacencyIterator begin, end;
	boost::tie(begin,end) = boost::adjacent_vertices(v,ag);
	for( AssemblyGraph::AdjacencyIterator z = begin; z != end; z++ )
	{
		this->initGraphDFS(ag, *z, v, colors, ag2cg);
	}
}


void
CompactAssemblyGraph::initGraph( const AssemblyGraph &ag )
{
    this->clear();
	this->_num_vertices = 0;

	size_t ag_vertices = boost::num_vertices(ag);

    //boost::dynamic_bitset<> colors( ag_vertices ); //std::vector<bool> colors( ag_vertices, 0 );
	//std::vector<Vertex> ag2cg( ag_vertices, 0 ); // associate each vertex of ag to a vertex of this CompactAssemblyGraph

	boost::dynamic_bitset<> *colors = new boost::dynamic_bitset<>(ag_vertices);
	std::vector<Vertex> *ag2cg = new std::vector<Vertex>(ag_vertices,0);

	AssemblyGraph::VertexIterator vbegin,vend;
	boost::tie(vbegin,vend) = boost::vertices(ag);

	for( AssemblyGraph::VertexIterator v = vbegin; v != vend; v++ )
	{
        if( boost::in_degree(*v,ag) == 0 && !colors->test(*v) )
        {
            Vertex new_v = boost::add_vertex(*this);
            std::list<Block> new_list = std::list<Block>(1,ag.getBlock(*v));
            _blockVector.push_back( new_list );

			colors->set(*v);
            ag2cg->at(*v) = new_v;

			_num_vertices++;

            AssemblyGraph::AdjacencyIterator begin, end;
            boost::tie(begin,end) = boost::adjacent_vertices(*v,ag);
            for( AssemblyGraph::AdjacencyIterator z = begin; z != end; z++ )
			{
				this->initGraphDFS( ag, *z, *v, colors, ag2cg );
			}
        }
	}

	delete colors;
	delete ag2cg;
}


void
CompactAssemblyGraph::computeEdgeWeights( MultiBamReader &masterBamReader, MultiBamReader &masterMpBamReader,
										  MultiBamReader &slaveBamReader, MultiBamReader &slaveMpBamReader )
{
	EdgeIterator ebegin,eend;
	boost::tie(ebegin,eend) = boost::edges(*this);

	for( EdgeIterator e = ebegin; e != eend; e++ )
	{
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), *this, *e);
		EdgeKindType kind = edge_prop.kind;

		std::list<Block>& b1 = _blockVector.at( boost::source(*e,*this) );
		std::list<Block>& b2 = _blockVector.at( boost::target(*e,*this) );

		switch(kind)
		{
			case MASTER_EDGE:
				this->getRegionScore(
					masterBamReader, masterMpBamReader, MASTER_EDGE, b1, b2,
					edge_prop.weight, edge_prop.rnum, edge_prop.min_cov
				);
				break;

			case SLAVE_EDGE:
				this->getRegionScore(
					slaveBamReader, slaveMpBamReader, SLAVE_EDGE, b1, b2,
					edge_prop.weight, edge_prop.rnum, edge_prop.min_cov
				);
				break;

			default:
				edge_prop.weight = 0.0;
				edge_prop.rnum = 0;
				edge_prop.min_cov = false;
				break;
		}

		// put edge weight
		boost::put( boost::edge_kind_t(), *this, *e, edge_prop );
	}
}


void CompactAssemblyGraph::getRegionScore( MultiBamReader &peBamReader, MultiBamReader &mpBamReader, EdgeKindType kind,
										   std::list<Block>& b1, std::list<Block>& b2,
										   double &weight, int32_t &rnum, bool &min_cov )
{
	std::vector< std::pair<double,int32_t> > mpStats, peStats;

	double mp_weight, pe_weight;
	int32_t mp_rnum, pe_rnum;
	bool mp_min_cov, pe_min_cov;

	if(peBamReader.size() > 0) getLibRegionScore( peBamReader, kind, b1, b2, pe_weight, pe_rnum, pe_min_cov );
	if(mpBamReader.size() > 0) getLibRegionScore( mpBamReader, kind, b1, b2, mp_weight, mp_rnum, mp_min_cov );

	min_cov = (pe_min_cov || mp_min_cov);

	// min number of evidences only for PE library
	if( pe_rnum >= 10 && mp_rnum < 10 ){ weight = pe_weight; rnum = pe_rnum; return; }
	// min number of evidences only for MP library
	if( mp_rnum >= 10 && pe_rnum < 10 ){ weight = mp_weight; rnum = mp_rnum; return; }
	// not enough evidences for both PE/MP libraries
	if( pe_rnum < 10 && mp_rnum < 10 ){ weight = -5.0; rnum = 0; return; }

	// enough evidences for both PE/MP libraries

	if( pe_weight >= 0 && mp_weight < 0 ){ weight = pe_weight; rnum = pe_rnum; return; }
	if( mp_weight >= 0 && pe_weight < 0 ){ weight = mp_weight; rnum = mp_rnum; return; }
	if( pe_weight < 0 && mp_weight < 0 ){ weight = -10.0; rnum = 0; return; }

	weight = pe_weight > mp_weight ? pe_weight : mp_weight;
	rnum = pe_weight > mp_weight ? pe_rnum : mp_rnum;

	return;
}


void CompactAssemblyGraph::getLibRegionScore( MultiBamReader &bamReader, EdgeKindType kind, std::list<Block>& b1, std::list<Block>& b2,
											  double &weight, int32_t &rnum, bool &min_cov )
{
	weight = -4;
	rnum = 0;
	min_cov = false;

	int32_t id, seq_len, start, end, region, s1, s2, t;
	uint64_t good_reads, exp_reads, num_reads;

	std::vector<double> score( bamReader.size(), -4 );
	std::vector<int32_t> r_num( bamReader.size(), 0 );
	std::vector<bool> cov( bamReader.size(), false );
	const RefVector& ref = bamReader.GetReferenceData();

	// this shouldn't happen
	if( kind != MASTER_EDGE && kind != SLAVE_EDGE ) return;
	if( b1.size() == 0 || b2.size() == 0 ) return;

	Frame& f1 = (kind == MASTER_EDGE) ? b1.front().getMasterFrame() : b1.front().getSlaveFrame();
	Frame& f2 = (kind == MASTER_EDGE) ? b2.front().getMasterFrame() : b2.front().getSlaveFrame();
	Frame& l1 = (kind == MASTER_EDGE) ? b1.back().getMasterFrame() : b1.back().getSlaveFrame();
	Frame& l2 = (kind == MASTER_EDGE) ? b2.back().getMasterFrame() : b2.back().getSlaveFrame();

	int32_t r1_beg = std::min( f1.getBegin(), l1.getBegin() );
	int32_t r1_end = std::max( f1.getEnd(), l1.getEnd() );
	int32_t r2_beg = std::min( f2.getBegin(), l2.getBegin() );
	int32_t r2_end = std::max( f2.getEnd(), l2.getEnd() );

	// skip included frames
	if( (r1_beg <= r2_beg && r1_end >= r2_end) ||
		(r2_beg <= r1_beg && r2_end >= r1_end) )
	{
		weight = -1;
		rnum = 0;
		min_cov = false;

		return;
	}

	int32_t gap = (r1_beg <= r2_beg) ? (r2_beg - r1_end + 1) : (r1_beg - r2_end + 1);

	// COMPUTE STATISTICS FOR EACH LIBRARY
	for( int lib=0; lib < bamReader.size(); lib++ )
	{
		int32_t isizeLibMean = bamReader.getISizeMean(lib);
		int32_t isizeLibStd = bamReader.getISizeStd(lib);
		int32_t coverageLib = bamReader.getCoverage(lib);

		int32_t minInsert = isizeLibMean - 3*isizeLibStd;
		int32_t maxInsert = isizeLibMean + 3*isizeLibStd;

		if(minInsert < 0) minInsert = 0;

		id = f1.getContigId();
		seq_len = ref[id].RefLength;

		t = (r1_beg <= r2_beg) ? (gap >= 0 ? r2_beg : r1_end) : (gap >= 0 ? r1_beg : r2_end);
		s1 = std::max( t - maxInsert, 0 );
		s2 = (r1_beg <= r2_beg) ? (gap >= 0 ? r1_end : r2_beg) : (gap >= 0 ? r2_end : r1_beg);

		if( seq_len - s1 < maxInsert )
		{
			weight = -2;
			continue;
		}

		if( gap >= maxInsert || s2 < s1 )
		{
			weight = -3;
			continue;
		}

		int32_t region = s2 - s1 + 1;
		std::vector<uint32_t> coverage( region, 0 );

		// retrieve BAM readers for current library
		bamReader.lockBamReader(lib);

		BamReader *reader = bamReader.getBamReader(lib);
		reader->SetRegion( id, s1, id, s2+1 );

		int32_t nh, xt;
		good_reads = 0;
		exp_reads = 0;
		num_reads = 0;

		BamAlignment align;
		while( reader->GetNextAlignment(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || align.Position < 0 || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;
			//if( !align.IsMateMapped() || align.RefID != align.MateRefID || align.MatePosition < t ) continue;

			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;

			for( int32_t i=startRead; i <= endRead; i++ ) if( i >= s1 && i <= s2 ) coverage[i-s1]++;

			if( !align.IsPaired() ) continue;

			//align.BuildCharData(); // fill string fields

			// if not defined, I assume read's multiplicity is 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard field
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';	// bwa field
			bool uniqMapRead = g_options.noMultiplicityFilter || (nh == 1 && xt == 'U');
			
			if( !uniqMapRead ) continue; // discard reads with multiplicity greater than 1

			int32_t startMate = align.MatePosition;
			int32_t endMate = startMate + readLength - 1;

			// don't count reads not completely included in the region
			if( startRead < s1 || startRead > s2 ) continue; //|| endRead > s2 ) continue;

			if( !align.IsReverseStrand() )
			{
				int32_t minInsertPos = startRead + minInsert;
				int32_t maxInsertPos = startRead + maxInsert;
				int32_t readOverlap = endRead > s2 ? s2-startRead+1 : readLength;

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

		// check coverage in the region
		for( size_t i=0; i < coverage.size(); i++ )
		{
			if( coverage[i] * 3 < coverageLib ) cov[lib] = false;
		}

		if( num_reads < 10 || exp_reads == 0 )
		{
			score[lib] = -5;
			r_num[lib] = 0;
		}
		else
		{
			score[lib] = good_reads / ((double)exp_reads);
			r_num[lib] = num_reads;
		}
	}

	// output statistics gained with the library with most evidences
	for( size_t i=0; i < score.size(); i++ )
	{
		if( i==0 )
		{
			weight = score[i];
			rnum = r_num[i];
			min_cov = cov[i];
		}
		else
		{
			if( r_num[i] > rnum ){ weight = score[i]; rnum = r_num[i]; }
			min_cov = min_cov || cov[i];
		}
	}
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

void CompactAssemblyGraph::bubbleDFS( Vertex v, std::vector<char> &colors, bool &found )
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
