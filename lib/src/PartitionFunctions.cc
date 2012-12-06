#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/filesystem.hpp>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "OptionsMerge.hpp"

#include "PartitionFunctions.hpp"
#include "graphs/PairedGraph.code.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "graphs/CompactAssemblyGraph.hpp"
#include "pctg/PairedContig.hpp"

#include "UtilityFunctions.hpp"

extern OptionsMerge g_options;

extern std::ofstream _g_statsFile;

extern MultiBamReader masterBam;
extern MultiBamReader masterMpBam;
extern MultiBamReader slaveBam;
extern MultiBamReader slaveMpBam;

std::list< std::vector<Block> >
partitionBlocks( const std::list<Block> &blocks )
{
    typedef AssemblyGraph::vertex_iterator VertexIterator;
    typedef AssemblyGraph::vertex_descriptor Vertex;
    typedef AssemblyGraph::adjacency_iterator AdjacencyIterator;
    typedef AssemblyGraph::inv_adjacency_iterator InvAdjacencyIterator;

    std::list< std::vector<Block> > blocksList;

    // group blocks between contigs which may be merged together through weaving
    std::vector< std::list<Block> > pairedContigsBlocks = partitionBlocksByPairedContigs( blocks );

    int z = 1; // assembly graph counter

    uint32_t ag_forks = 0, ag_linears = 0, ag_cycles = 0, ag_bubbles = 0; // counters for the different types of assemblies's graphs.

    // for each partition of blocks
    std::vector< std::list<Block> >::iterator pcb;
    for( pcb = pairedContigsBlocks.begin(); pcb != pairedContigsBlocks.end(); ++pcb )
    {
        // create an assembly graph
        AssemblyGraph ag( *pcb );
		//ag.computeEdgeWeights( masterBam, masterMpBam, slaveBam, slaveMpBam );
		
		//std::cerr << "AG " << z << std::endl;
		bool is_linear = false;
		
		while( !is_linear )
		{
			try
			{
				// check if graph still contains cycles (without forks).
				std::vector< size_t > ts;
				boost::topological_sort( ag, std::back_inserter(ts) );
				
				bool has_bubbles = ag.hasBubbles();
				
				std::vector<Vertex> v_forks;
				
				if( has_bubbles )
				{
					ag_bubbles++;
					
					std::stringstream ff1;
					ff1 << "./gam_graphs/AssemblyGraph_" << z << "_bubbles.dot";
					
					boost::filesystem::path p(ff1.str().c_str());
					if( not boost::filesystem::exists(p) )
					{
						std::ofstream ss1( ff1.str().c_str() );
						ag.writeGraphviz(ss1);
						ss1.close();
					}
				}

                ag.getForks( v_forks ); //if( !has_bubbles ) ag.getForks( v_forks );

                if( v_forks.size() == 0 ) is_linear = true;

                // Assembly graph with forks
                if( !has_bubbles && v_forks.size() >= 1 )
                {
                    ag_forks++;

                    std::stringstream ff1;
                    ff1 << "./gam_graphs/AssemblyGraph_" << z;
                    if( v_forks.size() > 1 ) ff1 << "_forks.dot"; else ff1 << "_fork.dot";

                    boost::filesystem::path p(ff1.str().c_str());
                    if( not boost::filesystem::exists(p) )
                    {
                        std::ofstream ss1( ff1.str().c_str() );
                        ag.writeGraphviz(ss1);
                        ss1.close();
                    }

                    std::stringstream ff2;
                    ff2 << "./gam_graphs/CompactGraph_" << z;
                    if( v_forks.size() > 1 ) ff2 << "_forks.dot"; else ff2 << "_fork.dot";

                    boost::filesystem::path p2(ff2.str().c_str());
                    if( not boost::filesystem::exists(p2) )
                    {
                        std::ofstream ss2( ff2.str().c_str() );
                        CompactAssemblyGraph cg(ag);
						//cg.computeEdgeWeights( masterBam, masterMpBam, slaveBam, slaveMpBam );
                        cg.writeGraphviz(ss2);
                        ss2.close();
                    }
                }

				// FINE GESTIONE MIS-ASSEMBLY/REPEATS

                if( !has_bubbles && v_forks.size() == 0 && boost::num_vertices(ag) > 0 )
                {
                    ag_linears++;

                    std::stringstream ff2;
                    ff2 << "./gam_graphs/AssemblyGraph_" << z << "_linear.dot";

                    boost::filesystem::path p(ff2.str().c_str());
                    if( not boost::filesystem::exists(p) )
                    {
                        std::ofstream ss2( ff2.str().c_str() );
                        ag.writeGraphviz(ss2);
                        ss2.close();
                    }

                    std::stringstream ff3;
                    ff3 << "./gam_graphs/CompactGraph_" << z << "_linear.dot";

                    boost::filesystem::path p2(ff3.str().c_str());
                    if( not boost::filesystem::exists(p2) )
                    {
                        std::ofstream ss3( ff3.str().c_str() );
                        CompactAssemblyGraph cg(ag);
						//cg.computeEdgeWeights( masterBam, masterMpBam, slaveBam, slaveMpBam );
                        cg.writeGraphviz(ss3);
                        ss3.close();
                    }

                    blocksList.push_back(ag.getBlocksVector());
                }
                else if( !has_bubbles && v_forks.size() >= 1 && boost::num_vertices(ag) > 0 )
                {
                    blocksList.push_back(ag.getBlocksVector());
                    break;
                }
                else
                {
					// grafi con bubbles sono eventualmente scartati dopo
					blocksList.push_back(ag.getBlocksVector());
					
                    break; //WARNING => temporaneo: solo per fare uscire dal while senza risolvere bubbles
                    ag.removeForks();
                }
            }
            catch( boost::not_a_dag ) // if the graph is not a DAG, remove cycles.
            {
                ag_cycles++;

                /* DEBUG - Assembly graph with cycles */
                std::stringstream ff1;
                ff1 << "./gam_graphs/AssemblyGraph_" << z << "_cyclic.dot";

                boost::filesystem::path p(ff1.str().c_str());
                if( not boost::filesystem::exists(p) )
                {
                    std::ofstream ss1( ff1.str().c_str() );
                    ag.writeGraphviz(ss1);
                    ss1.close();
                }
                /* END DEBUG*/

                break; //WARNING => temporaneo: solo per fare uscire dal while senza risolvere i cicli

                ag.removeCycles();
                //cyclesRemoved = true;
            }
        } // end of while( !is_linear )

        z++; // increase assembly graph counter
    }

    _g_statsFile << "[graphs stats]\n"
		<< "Linears = " << ag_linears << "\n"
		<< "Forks = " << ag_forks << "\n"
		<< "Bubbles = " << ag_bubbles << "\n"
		<< "Cyclics = " << ag_cycles << "\n"
		<< std::endl;

    return blocksList;
}


std::vector<double> computeZScore( MultiBamReader &multiBamReader, const uint64_t &refID, uint32_t start, uint32_t end )
{
	uint32_t libs;
	BamReader* bamReader;
	double lib_isize_mean, lib_isize_std;

	uint32_t minInsertNum = 5;

	libs = multiBamReader.size();
	std::vector<double> z_score( libs, 0.0 );

	if( libs == 0 ) return z_score;

	unsigned int times_std = 3;

	for( int i=0; i<libs; i++ )
	{
		lib_isize_mean = multiBamReader.getISizeMean(i);
		lib_isize_std = multiBamReader.getISizeStd(i);
		bamReader = multiBamReader.getBamReader(i);

		if( lib_isize_std == 0 ) continue;

		uint32_t min_insert = ( lib_isize_mean > times_std*lib_isize_std ) ? lib_isize_mean - times_std*lib_isize_std : 0;
		uint32_t max_insert = lib_isize_mean + times_std*lib_isize_std;

		BamAlignment align;
		uint64_t inserts=0, spanCov=0;
		int32_t nh, xt; // molteplicità delle read (nh->standard, xt->bwa)

		multiBamReader.lockBamReader(i);

		bamReader->SetRegion( refID, start, refID, end+1 );
		while( bamReader->GetNextAlignmentCore(align) ) // for each read in the region
		{
			if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ||
				!align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

			int32_t read_start = align.Position;
			int32_t read_end = align.GetEndPosition() - 1;
			int32_t read_len = read_end - read_start + 1;
			int32_t mate_start = align.MatePosition;
			int32_t mate_end = align.MatePosition + read_len - 1;

			if( read_start < start || read_end > end ) continue;
			if( mate_start < start || mate_end > end ) continue;

			align.BuildCharData();
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard SAM format field
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // bwa field
			bool is_uniq_mapped = (nh == 1 && xt == 'U');

			if( !is_uniq_mapped ) continue;

			if( align.IsFirstMate() )
			{
				if( read_start < mate_start )
				{
					int32_t i_size = (mate_start + read_len) - read_start;
					if( i_size < min_insert || i_size > max_insert ) continue;

					inserts++;
					spanCov += i_size;
				}
				else
				{
					int32_t i_size = read_end - mate_start + 1;
					if( i_size < min_insert || i_size > max_insert ) continue;

					inserts++;
					spanCov += i_size;
				}
			}
		} // end while

		multiBamReader.unlockBamReader(i);

		if( inserts > minInsertNum )
		{
			double localMean = spanCov/(double)inserts;
			z_score[i] = (localMean - lib_isize_mean)/(double)(lib_isize_std/sqrt(inserts));
		}
	}

	return z_score;
}


std::vector< std::list<Block> >
partitionBlocksByPairedContigs( const std::list< Block > &blocks )
{
    typedef PairedContigGraph<> PCGraph;
    typedef boost::graph_traits< PCGraph >::vertex_descriptor Vertex;

    PCGraph pcg(blocks);

    // compute connected components
    std::vector< uint64_t > component( boost::num_vertices(pcg) );
    int num = boost::connected_components( pcg, &component[0] );

    // each connected component consists of blocks of contigs which may be extended through weaving.
    std::vector< std::list< Block > > pairedContigs(num);
    for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); ++b )
    {
        Vertex v = pcg.getMasterVertex(*b);
        pairedContigs.at( component.at(v) ).push_back(*b);
    }

    // print pairedcontigs graphs
    for( UIntType i=0; i < pairedContigs.size(); i++ )
    {
        PCGraph g( pairedContigs[i] );

        std::stringstream ss;
        ss << "./gam_graphs/ContigGraph_" << i << ".dot";
        std::ofstream ofs( ss.str().c_str() );

        g.writeGraphviz(ofs);

        ofs.close();
    }

    return pairedContigs;
}
