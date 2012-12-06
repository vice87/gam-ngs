/*
 * File:   AssemblyGraph.cc
 * Author: vice
 *
 * Created on 25 maggio 2011, 12.07
 */

#include "graphs/AssemblyGraph.hpp"

#include <vector>
#include <iostream>
#include <iomanip>

#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>

#include "strand_fixer/RelativeStrand.hpp"
#include "OrderingFunctions.hpp"


AssemblyGraph::AssemblyGraph()
{
    std::list< Block > blocks;
    this->initGraph( blocks );
}


AssemblyGraph::AssemblyGraph( const std::list< Block > &blocks )
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


/*std::list<Block>
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
}*/


void //std::vector< std::vector<Block> >
AssemblyGraph::removeCycles()
{
	std::vector< uint64_t > component( boost::num_vertices(*this) );
	std::vector< uint64_t > discoverTime( boost::num_vertices(*this) );
	std::vector< uint64_t > root( boost::num_vertices(*this) );
	std::vector< boost::default_color_type > color( boost::num_vertices(*this) );

    // compute strongly connected components
    uint64_t sccNum = boost::strong_components( *this, &component[0],
            boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discoverTime[0]) );

    // group vertices by scc
    std::vector< std::list<uint64_t> > components(sccNum);
    for( uint64_t h=0; h < component.size(); ++h )
        components.at( component.at(h) ).push_back(h);

    // compute the number of strong components which include more than one vertex
    UIntType scc = 0;
    for( uint64_t h=0; h != sccNum; ++h ) if( components[h].size() > 1 ) scc++;

    std::vector< std::list<Block> > sComponents(scc);
    scc = 0;

    // fill strong components data structure
    for( uint64_t h=0; h < sccNum; ++h )
    {
        if( components[h].size() > 1 )
        {
			for( std::list<uint64_t>::iterator j = components[h].begin(); j != components[h].end(); j++ )
                sComponents[scc].push_back( this->_blockVector.at(*j) );

            scc++;
        }
    }

    // remove strong components from graph
    std::list<Block> newBlockList;
    for( UIntType i = 0; i < this->_blockVector.size(); i++ )
    {
        if( components.at( component.at(i) ).size() == 1 )
			newBlockList.push_back( this->_blockVector[i] );
    }

    // re-initialize the graph with blocks that do not create cycles
    this->initGraph( newBlockList );

    //std::cout << "Filtered blocks after cycles removal = " << newBlockVector.size() << std::endl << std::flush;

    //return sComponents;
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


void
AssemblyGraph::computeEdgeWeights( MultiBamReader &masterBamReader, MultiBamReader &masterMpBamReader, MultiBamReader &slaveBamReader, MultiBamReader &slaveMpBamReader )
{
	EdgeIterator ebegin,eend;
	boost::tie(ebegin,eend) = boost::edges(*this);

	for( EdgeIterator e = ebegin; e != eend; e++ )
	{
		//EdgeKindType kind = boost::get(boost::edge_kind_t(), *this, *e);
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), *this, *e);
		EdgeKindType kind = edge_prop.kind;

		Block& b1 = _blockVector.at( boost::source(*e,*this) );
		Block& b2 = _blockVector.at( boost::target(*e,*this) );

		switch(kind)
		{
			case MASTER_EDGE:
				edge_prop.weight = this->getRegionScore( masterBamReader, masterMpBamReader, MASTER_EDGE, b1, b2 );
				break;

			case SLAVE_EDGE:
				edge_prop.weight = this->getRegionScore( slaveBamReader, slaveMpBamReader, SLAVE_EDGE, b1, b2 );
				break;

			case BOTH_EDGE:
				edge_prop.weight = 0.0;
				break;
		}

		// put edge weight
		boost::put( boost::edge_kind_t(), *this, *e, edge_prop );
	}
}


double AssemblyGraph::getRegionScore( MultiBamReader &peBamReader, MultiBamReader &mpBamReader, EdgeKindType kind, Block& b1, Block& b2 )
{
	//std::cerr << "PE" << std::endl;
	//std::vector<double> peScore = getLibRegionScore( peBamReader, kind, b1, b2 );
	//std::cerr << "MP" << std::endl;
	std::vector<double> mpScore = getLibRegionScore2( mpBamReader, kind, b1, b2 );
	
	return mpScore[0];

	//double maxScore = 0.0;
	//for( size_t i=0; i < peScore.size(); i++ ) if( peScore[i] > maxScore ) maxScore = peScore[i];
	//for( size_t i=0; i < mpScore.size(); i++ ) if( mpScore[i] > maxScore ) maxScore = mpScore[i];

	//return maxScore;
}


std::vector<double> AssemblyGraph::getLibRegionScore2( MultiBamReader &bamReader, EdgeKindType kind, Block& b1, Block& b2 )
{
	int32_t id, start, end, gap, region, s1, s2, t;
	uint64_t mates;
	
	std::vector<double> score( bamReader.size(), -4 );
	const RefVector& ref = bamReader.GetReferenceData();
	
	// this shouldn't happen
	if( kind != MASTER_EDGE && kind != SLAVE_EDGE ) return score;
	
	// COMPUTE STATISTICS FOR EACH LIBRARY
	for( int lib=0; lib < bamReader.size(); lib++ )
	{
		int32_t isizeLibMean = bamReader.getISizeMean(lib);
		int32_t isizeLibStd = bamReader.getISizeStd(lib);
		
		int32_t minInsert = isizeLibMean - 3*isizeLibStd;
		int32_t maxInsert = isizeLibMean + 3*isizeLibStd;
		
		if(minInsert < 0) minInsert = 0;
		
		Frame& f1 = (kind == MASTER_EDGE) ? b1.getMasterFrame() : b1.getSlaveFrame();
		Frame& f2 = (kind == MASTER_EDGE) ? b2.getMasterFrame() : b2.getSlaveFrame();
		
		id = f1.getContigId();
		
		start = (f1.getBegin() <= f2.getBegin()) ? std::max( f1.getEnd() - maxInsert, 0 ) : std::max( f2.getEnd() - maxInsert, 0 ); //std::min( f1.getBegin(), f2.getBegin() );
		end = std::max( f1.getEnd(), f2.getEnd() );
		region = end - start + 1;
		
		// skip overlapping frames
		if( (f1.getBegin() <= f2.getBegin() && f1.getEnd() >= f2.getBegin()) || 
			(f2.getBegin() <= f1.getBegin() && f2.getEnd() >= f1.getBegin()) )
		{
			score[lib] = -1;
			continue;
		}
		
		gap = (f1.getBegin() <= f2.getBegin()) ? (f2.getBegin() - f1.getEnd() + 1) : (f1.getBegin() - f2.getEnd() + 1);
		
		if( region < maxInsert )
		{
			score[lib] = -2;
			continue;
		}
		
		if( gap >= maxInsert )
		{
			score[lib] = -3;
			continue;
		}
		
		s1 = start; //(f1.getBegin() <= f2.getBegin()) ? std::max( f1.getEnd() - maxInsert, 0 ) : std::max( f2.getEnd() - maxInsert, 0 );
		s2 = (f1.getBegin() <= f2.getBegin()) ? f1.getEnd() : f2.getEnd();
		t = (f1.getBegin() <= f2.getBegin()) ? f2.getBegin() : f1.getBegin();
		
		// retrieve BAM readers for current library
		bamReader.lockBamReader(lib);
		
		BamReader *reader = bamReader.getBamReader(lib);
		reader->SetRegion( id, s1, id, s2+1 );
		
		int32_t nh, xt;
		mates = 0;
		
		BamAlignment align;
		while( reader->GetNextAlignmentCore(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;
			if( !align.IsMateMapped() || align.RefID != align.MateRefID || align.MatePosition < t ) continue;
			
			align.BuildCharData(); // fill string fields
			
			// if not defined, I assume read's multiplicity is 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard field
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';	// bwa field
			if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1
			
			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;
			int32_t startMate = align.MatePosition;
			
			// don't count reads not completely included in the region
			if( startRead < s1 || endRead > s2 ) continue;
			
			if( align.IsMateMapped() && align.RefID == align.MateRefID && startMate >= t ) mates++;
		} // end while
		
		bamReader.unlockBamReader(lib);
		
		score[lib] = mates;
	}
	
	return score;
}


std::vector<double> AssemblyGraph::getLibRegionScore( MultiBamReader &bamReader, EdgeKindType kind, Block& b1, Block& b2 )
{
	int32_t id, start, end, seq_len, ltail, rtail;

	std::vector<double> score( bamReader.size(), -4 );
	const RefVector& ref = bamReader.GetReferenceData();
	
	if( kind != MASTER_EDGE && kind != SLAVE_EDGE ) return score;

	// COMPUTE STATISTICS FOR EACH LIBRARY
	for( int lib=0; lib < bamReader.size(); lib++ )
	{
		int32_t isizeLibMean = bamReader.getISizeMean(lib);
		int32_t isizeLibStd = bamReader.getISizeStd(lib);

		int32_t minInsert = isizeLibMean - 3*isizeLibStd;
		int32_t maxInsert = isizeLibMean + 3*isizeLibStd;
		
		Frame& f1 = (kind == MASTER_EDGE) ? b1.getMasterFrame() : b1.getSlaveFrame();
		Frame& f2 = (kind == MASTER_EDGE) ? b2.getMasterFrame() : b2.getSlaveFrame();

		// one frame included into the other one
		if( (f1.getBegin() <= f2.getBegin() && f2.getBegin() <= f1.getEnd()) || (f2.getBegin() <= f1.getBegin() && f1.getBegin() <= f2.getEnd()) )
		{
			score[lib] = -1;
			continue;
		}
		
		id = f1.getContigId();
		seq_len = ref[id].RefLength;
		
		if( f1.getBegin() <= f2.getBegin() )
		{
			rtail = seq_len - f2.getBegin();
			ltail = f1.getEnd() + 1;
			
			start = (rtail >= ltail) ? std::max( f1.getBegin(), f1.getEnd() - isizeLibMean ) : f1.getEnd();
			end = (rtail >= ltail) ? f2.getBegin() : std::min( f2.getEnd(), f2.getBegin() + isizeLibMean );
		}
		else
		{
			rtail = seq_len - f1.getBegin();
			ltail = f2.getEnd() + 1;
			
			start = (rtail >= ltail) ? std::max( f2.getBegin(), f2.getEnd() - isizeLibMean ) : f2.getEnd();
			end = (rtail >= ltail) ? f1.getBegin() : std::min( f1.getEnd(), f1.getBegin() + isizeLibMean );
		}

		int32_t region = end - start + 1;
		if( region < 100 ){ score[lib] = -2; continue; } //if( isizeLibMean + 3*isizeLibStd > region ) continue;

		std::vector<uint64_t> allCoverage( region, 0 );
		std::vector<uint64_t> suspCoverage( region, 0 );

		// retrieve BAM readers for current library
		BamReader *reader = bamReader.getBamReader(lib);

		bamReader.lockBamReader(lib);
		
		if( rtail >= ltail ) this->updateRightRegionCoverage( bamReader.getBamReader(lib), id, start, end, isizeLibMean, isizeLibStd, allCoverage, suspCoverage );
			else this->updateLeftRegionCoverage( bamReader.getBamReader(lib), id, start, end, isizeLibMean, isizeLibStd, allCoverage, suspCoverage );
		
		bamReader.unlockBamReader(lib);

		// find putative suspicious zone
		uint64_t totAllCoverage = 0;
		uint64_t totSuspCoverage = 0;

		double meanAllCoverage;
		double meanSuspCoverage;

		double maxScore;

		uint32_t windowSize = isizeLibMean;
		uint32_t windowStep = 100;

		if(region < windowSize) // if region shorter than window size, only one window
		{
			for( size_t i=0; i < region; i++ )
			{
				totAllCoverage += allCoverage[i];
				totSuspCoverage += suspCoverage[i];
			}

			meanAllCoverage = totAllCoverage / (double)region; // this is the "window" total coverage
			meanSuspCoverage = totSuspCoverage/ (double)region; // this is the "window" suspicious coverage

			maxScore = ( meanAllCoverage != 0 ) ? meanSuspCoverage / meanAllCoverage : 0.0;
			//maxScore = ( meanSuspCoverage >= 0 && meanAllCoverage > 0 ) ? meanSuspCoverage / meanAllCoverage : 0.0;
		}
		else //otherwise compute score on sliding window of 200 bp
		{
			unsigned int startWindow = 0;
			unsigned int endWindow = windowSize;
			unsigned int winSize   = windowSize;
			bool last = false;

			// first window
			for(size_t i = startWindow; i < endWindow; i++ )
			{
				totAllCoverage += allCoverage[i];
				totSuspCoverage += suspCoverage[i];
			}

			meanAllCoverage = totAllCoverage/(double)winSize;
			meanSuspCoverage = totSuspCoverage/(double)winSize;

			maxScore = ( meanAllCoverage != 0 ) ? meanSuspCoverage / meanAllCoverage : 0.0;
			//maxScore = ( meanSuspCoverage >= 0 && meanAllCoverage > 0 ) ? meanSuspCoverage / meanAllCoverage : 0.0;

			//now update
			startWindow += windowStep;
			endWindow += windowStep;

			if( endWindow > region ) endWindow = region;

			// inner windows
			while( endWindow < region )
			{
				totAllCoverage = 0;
				totSuspCoverage = 0;

				for(size_t i = startWindow; i < endWindow; i++ )
				{
					totAllCoverage += allCoverage[i];
					totSuspCoverage += suspCoverage[i];
				}

				meanAllCoverage = totAllCoverage/(double)(endWindow - startWindow);
				meanSuspCoverage = totSuspCoverage/(double)(endWindow - startWindow);

				maxScore = ( meanAllCoverage != 0 ) ? meanSuspCoverage / meanAllCoverage : 0.0;
				//maxScore = std::max( maxScore, (meanSuspCoverage >= 0 && meanAllCoverage > 0) ? meanSuspCoverage / meanAllCoverage : 0.0 );

				startWindow += windowStep;
				endWindow += windowStep;
				if( endWindow > region ) endWindow = region;
			}

			// last window
			totAllCoverage = 0;
			totSuspCoverage = 0;

			for(size_t i = startWindow; i < endWindow; i++ )
			{
				totAllCoverage += allCoverage[i];
				totSuspCoverage += suspCoverage[i];
			}

			meanAllCoverage = totAllCoverage/(double)(endWindow - startWindow);
			meanSuspCoverage = totSuspCoverage/(double)(endWindow - startWindow);

			maxScore = ( meanAllCoverage != 0 ) ? meanSuspCoverage / meanAllCoverage : 0.0;
			//maxScore = std::max( maxScore, (meanSuspCoverage >= 0 && meanAllCoverage > 0) ? meanSuspCoverage / meanAllCoverage : 0.0 );
		}

		score[lib] = maxScore;
	}

	return score;
}


void AssemblyGraph::updateRightRegionCoverage( BamReader *bamReader, int32_t id, int32_t start, int32_t end,
										  double isizeLibMean, double isizeLibStd,
										  std::vector<uint64_t> &allCoverage, std::vector<uint64_t> &suspCoverage )
{
	bamReader->SetRegion( id, start, id, end+1 );

	int32_t nh, xt;
	int32_t ctgLen = (bamReader->GetReferenceData()).at(id).RefLength;

	BamAlignment align;
	while( bamReader->GetNextAlignmentCore(align) )
	{
		// discard bad quality reads
		if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

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

		int32_t maxInsert = isizeLibMean + 3*isizeLibStd;
		int32_t minInsert = isizeLibMean - 3*isizeLibStd;
		
		if(minInsert < 0) minInsert = 0;

		// don't count reads not completely included in the region
		if( startRead < start || endRead > end ) continue;

		// reads with mate mapped in the same sequence
		if( align.IsMateMapped() && align.RefID == align.MateRefID )
		{
			//if( (align.IsFirstMate() && !align.IsReverseStrand()) || (align.IsSecondMate() && align.IsReverseStrand()) )
			if( !align.IsReverseStrand() )
			{
				//int32_t endPos = startRead + minInsert;
				//int32_t endPos2 = startRead + maxInsert;
				
				if( align.IsMateReverseStrand() ) // correct orientation of the pair
				{
					for( int32_t i=startRead; i <= endRead; i++ ) allCoverage[i-start]++;
				}
				else
				{
					for( int32_t i=startRead; i <= endRead; i++ )
					{
						suspCoverage[i-start]++;
						allCoverage[i-start]++;
					}
				}
			}
		}

		// mate/pair unmapped or mapped to different sequence
		if( !align.IsMateMapped() || (align.IsMateMapped() && align.RefID != align.MateRefID) )
		{
			//if( (align.IsFirstMate() && !align.IsReverseStrand()) || (align.IsSecondMate() && align.IsReverseStrand()) )
			if( !align.IsReverseStrand() )
			{
				int32_t endPos = startRead + minInsert;
				int32_t endPos2 = startRead + maxInsert;

				if( (endPos < ctgLen && endPos2 < ctgLen) || !align.IsMateMapped() )
				{
					for( int32_t i=startRead; i <= endRead; i++ )
					{
						suspCoverage[i-start]++;
						allCoverage[i-start]++;
					}
				}
			}
		}
	} // end while
}


void AssemblyGraph::updateLeftRegionCoverage( BamReader *bamReader, int32_t id, int32_t start, int32_t end,
										  double isizeLibMean, double isizeLibStd,
										  std::vector<uint64_t> &allCoverage, std::vector<uint64_t> &suspCoverage )
{
	bamReader->SetRegion( id, start, id, end+1 );
	
	int32_t nh, xt;
	int32_t ctgLen = (bamReader->GetReferenceData()).at(id).RefLength;
	
	BamAlignment align;
	while( bamReader->GetNextAlignmentCore(align) )
	{
		// discard bad quality reads
		if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;
		
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
		
		int32_t maxInsert = isizeLibMean + 3*isizeLibStd;
		int32_t minInsert = isizeLibMean - 3*isizeLibStd;
		
		if(minInsert < 0) minInsert = 0;
		
		// don't count reads not completely included in the region
		if( startRead < start || endRead > end ) continue;
		
		// reads with mate mapped in the same sequence
		if( align.IsMateMapped() && align.RefID == align.MateRefID )
		{
			//if( (align.IsFirstMate() && !align.IsReverseStrand()) || (align.IsSecondMate() && align.IsReverseStrand()) )
			if( align.IsReverseStrand() )
			{
				//int32_t startPos = endRead - minInsert;
				//int32_t startPos2 = endRead - maxInsert;
				
				if( !align.IsMateReverseStrand() ) // correct orientation of the pair
				{
					for( int32_t i=startRead; i <= endRead; i++ )
						allCoverage[i-start]++;
				}
				else
				{
					for( int32_t i=startRead; i <= endRead; i++ )
					{
						allCoverage[i-start]++;
						suspCoverage[i-start]++;
					}
				}
				
			}
		}
		
		// mate/pair unmapped or mapped to different sequence
		if( !align.IsMateMapped() || (align.IsMateMapped() && align.RefID != align.MateRefID) )
		{
			//if( (align.IsFirstMate() && !align.IsReverseStrand()) || (align.IsSecondMate() && align.IsReverseStrand()) )
			if( align.IsReverseStrand() )
			{
				int32_t startPos = endRead - minInsert;
				int32_t startPos2 = endRead - maxInsert;
				
				if( (startPos >= 0 && startPos2 >= 0) || !align.IsMateMapped() )
				{
					for( int32_t i=startRead; i <= endRead; i++ )
					{
						suspCoverage[i-start]++;
						allCoverage[i-start]++;
					}
				}
			}
		}
	} // end while
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

void AssemblyGraph::getForks( std::vector<Vertex> &v_forks )
{
	typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
	typedef boost::graph_traits<AssemblyGraph>::edge_iterator EdgeIterator;

	VertexIterator vbegin,vend;
	boost::tie(vbegin,vend) = boost::vertices(*this);

	v_forks.clear();

	for (VertexIterator v=vbegin; v!=vend; v++)
	{
		int in_deg = boost::in_degree(*v,*this);
		int out_deg = boost::out_degree(*v,*this);
		if( in_deg > 1 ) v_forks.push_back(*v);
		if( out_deg > 1 ) v_forks.push_back(*v);
	}
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

bool AssemblyGraph::bubbleDFS( Vertex v, std::vector<char> &colors, bool &found )
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