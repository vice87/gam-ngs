#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/topological_sort.hpp>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "PartitionFunctions.hpp"
#include "graphs/PairedGraph.code.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "pctg/PairedContig.hpp"

extern std::ofstream _g_statsFile;

std::list< std::vector<Block> >
partitionBlocks( const std::vector<Block> &blocks, const Options &options )
{
    typedef AssemblyGraph::vertex_iterator VertexIterator;
    typedef AssemblyGraph::vertex_descriptor Vertex;
    typedef AssemblyGraph::adjacency_iterator AdjacencyIterator;
    typedef AssemblyGraph::inv_adjacency_iterator InvAdjacencyIterator;

    std::list< std::vector<Block> > blocksList;

    // group blocks between contigs which may be merged together through weaving
    std::vector< std::vector<Block> > pairedContigsBlocks =
            partitionBlocksByPairedContigs( blocks );

    int z = 1; // assembly graph counter

    uint32_t ag_forks = 0, ag_linears = 0, ag_cycles = 0; // counters for the different types of assemblies's graphs.

    // for each partition of blocks
    std::vector< std::vector<Block> >::iterator pcb;
    for( pcb = pairedContigsBlocks.begin(); pcb != pairedContigsBlocks.end(); pcb++ )
    {
        // create an assembly graph
        AssemblyGraph ag( *pcb );

        bool is_acyclic = false;

        try
        {
            // check if graph contains cycles.
            std::vector< size_t > ts;
            boost::topological_sort( ag, std::back_inserter(ts) );
            is_acyclic = true;

            // if the graph is acyclic, check how many forks it contains (nodes with in/out degree = 2)
            VertexIterator vbegin,vend;
            boost::tie(vbegin,vend) = boost::vertices(ag);
            uint32_t num_forks = 0;
            Vertex v_fork;
            bool fork_in = false;
            bool fork_out = false;

            for (VertexIterator v=vbegin; v!=vend; v++)
            {
                int in = boost::in_degree(*v,ag);
                int out = boost::out_degree(*v,ag);
                if( in > 1 ){ num_forks++; v_fork = *v; fork_in = true; }
                if( out > 1 ){ num_forks++; v_fork = *v; fork_out = true; }
            }

            // DEBUG - Assembly graph with forks
            if( num_forks >= 1 )
			{
				ag_forks++;

				std::stringstream ff1;
				ff1 << "./gam_graphs/AssemblyGraph_" << z;
				if( num_forks > 1 ) ff1 << "_forks.dot"; else ff1 << "_fork.dot";
				std::ofstream ss1( ff1.str().c_str() );
				ag.writeGraphviz(ss1);
				ss1.close();
			}
            // END of DEBUG

            // print mates information of the forking to discover mis-assemblies or repeats.
            /*if( num_forks == 1 )
            {
                AdjacencyIterator begin,end;
                InvAdjacencyIterator ibegin, iend;
                boost::tie(begin,end) = boost::adjacent_vertices(v_fork,ag);
                boost::tie(ibegin,iend) = boost::inv_adjacent_vertices(v_fork,ag);

                Block shared = ag.getBlock(v_fork);
                Block master, slave;

                if( fork_in )
                {
                    for( InvAdjacencyIterator v = ibegin; v != iend; v++ )
                    {
                        if( Block::shareMasterContig(shared,ag.getBlock(*v)) ) master = ag.getBlock(*v); else slave = ag.getBlock(*v);
                    }
                }
                else
                {
                    for( AdjacencyIterator v = begin; v != end; v++ )
                    {
                        if( Block::shareMasterContig(shared,ag.getBlock(*v)) ) master = ag.getBlock(*v); else slave = ag.getBlock(*v);
                    }
                }


                BamReader bamMaster,bamSlave;
                bamMaster.Open( options.masterBamFile );
                bamMaster.OpenIndex( options.masterBamFile + ".bai" );
                bamSlave.Open( options.slaveBamFiles.at( shared.getSlaveFrame().getAssemblyId() ) );
                bamSlave.OpenIndex( options.slaveBamFiles.at( shared.getSlaveFrame().getAssemblyId() ) + ".bai" );

                int32_t nh, xt;
                uint64_t m_pairs = 0, s_pairs = 0, m_mates = 0, s_mates = 0;

                // master pairs/mates information
                std::cerr << "Mates information for assembly graph " << z << std::endl;

                Frame m_frame1 = shared.getMasterFrame();
                Frame m_frame2 = master.getMasterFrame();
				Frame s_frame1 = shared.getSlaveFrame();
				Frame s_frame2 = slave.getSlaveFrame();

                int32_t id,s1,e1,s2,e2,len1,len2;
				int32_t sr,er,mr;

                // master interval computation
				if( m_frame1.getBegin() <= m_frame2.getBegin() )
				{
					len1 = bamMaster.GetReferenceData().at( m_frame2.getContigId() ).RefLength;
					s1 = m_frame1.getEnd();
					e1 = (m_frame1.getEnd() + 4000 < len1) ? m_frame1.getEnd() + 4000 : len1;

					if( m_frame2.getEnd() <= m_frame1.getEnd() ) std::cerr << "warning: one master frame is included into the other." << std::endl;
				}
				else
				{
					s1 = (m_frame1.getEnd() >= 4000) ? m_frame1.getEnd() - 4000 : 0;
					e1 = m_frame1.getEnd();

					if( m_frame1.getEnd() <= m_frame2.getEnd() ) std::cerr << "warning: one master frame is included into the other." << std::endl;
				}

				// slave interval computation
				if( s_frame1.getBegin() <= s_frame2.getBegin() )
				{
					len2 = bamSlave.GetReferenceData().at( s_frame2.getContigId() ).RefLength;
					s2 = s_frame1.getEnd();
					e2 = (s_frame1.getEnd() + 4000 < len2) ? s_frame1.getEnd() + 4000 : len2;

					if( s_frame2.getEnd() <= s_frame1.getEnd() ) std::cerr << "warning: one slave frame is included into the other." << std::endl;
				}
				else
				{
					s2 = (s_frame1.getEnd() >= 4000) ? s_frame1.getEnd() - 4000 : 0;
					e2 = s_frame1.getEnd();

					if( s_frame1.getEnd() <= s_frame2.getEnd() ) std::cerr << "warning: one slave frame is included into the other." << std::endl;
				}

				if( e2-s2 < e1-s1 )
				{
					if( m_frame1.getBegin() <= m_frame2.getBegin() ) e1 -= e1-s1-e2+s2; else s1 += e1-s1-e2+s2;
				}

				if( e2-s2 > e1-s1 )
				{
					if( s_frame1.getBegin() <= s_frame2.getBegin() ) e2 -= e2-s2-e1+s1; else s2 += e2-s2-e1+s1;
				}

                BamAlignment align;

				bamMaster.SetRegion( m_frame1.getContigId(), s1, m_frame1.getContigId(), e1 ); //bamMaster.SetRegion( id, s1, id, e2 );
                while( bamMaster.GetNextAlignmentCore(align) )
                {
                    if( !align.IsPaired() ) continue;
                    if( !align.IsMapped() || !align.IsMateMapped() ) continue;

                    if( align.RefID != align.MateRefID ) continue; // pairs must align in the same contig

                    align.BuildCharData(); // fill string fields

                    // se la molteplicità non è stata definita, assumo che sia pari ad 1
                    if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
                    if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';

                    if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

					if( m_frame1.getBegin() <= m_frame2.getBegin() && s1 <= align.Position && align.Position <= e1 && align.MatePosition < s1 ) m_pairs++;
					if( m_frame1.getBegin() > m_frame2.getBegin() && s1 <= align.Position && align.Position <= e1 && align.MatePosition > e1 ) m_pairs++;
					//if( sr <= align.Position && align.GetEndPosition()-1 <= er && align.MatePosition >= sr && align.MatePosition <= er ) m_pairs++;
                    //if( s1 <= align.Position && align.GetEndPosition()-1 <= e1 && align.MatePosition >= s2 && align.MatePosition <= e2 ) m_pairs++;
                    //if( s2 <= align.Position && align.GetEndPosition()-1 <= e2 && align.MatePosition >= s1 && align.MatePosition <= e1 ) m_pairs++;
                }

                // open mates alignment, if specified.
                if( options.masterMpBamFile != "" )
				{
					bamMaster.Open( options.masterMpBamFile );
					bamMaster.OpenIndex( options.masterMpBamFile + ".bai");

					bamMaster.SetRegion( m_frame1.getContigId(), s1, m_frame1.getContigId(), e1 ); //bamMaster.SetRegion( id, s1, id, e2 );
					while( bamMaster.GetNextAlignmentCore(align) )
					{
						if( !align.IsPaired() ) continue;
						if( !align.IsMapped() || !align.IsMateMapped() ) continue;

						if( align.RefID != align.MateRefID ) continue; // pairs must align in the same contig

						align.BuildCharData(); // fill string fields

						// se la molteplicità non è stata definita, assumo che sia pari ad 1
						if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
						if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';

						if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

						if( m_frame1.getBegin() <= m_frame2.getBegin() && s1 <= align.Position && align.Position <= e1 && align.MatePosition < s1 ) m_mates++;
						if( m_frame1.getBegin() > m_frame2.getBegin() && s1 <= align.Position && align.Position <= e1 && align.MatePosition > e1 ) m_mates++;
						//if( sr <= align.Position && align.GetEndPosition()-1 <= er && align.MatePosition >= sr && align.MatePosition <= er ) m_mates++;
						//if( s1 <= align.Position && align.GetEndPosition()-1 <= e1 && align.MatePosition >= s2 && align.MatePosition <= e2 ) m_mates++;
						//if( s2 <= align.Position && align.GetEndPosition()-1 <= e2 && align.MatePosition >= s1 && align.MatePosition <= e1 ) m_mates++;
					}
				} // end of mates number computation

				std::cerr << "master interval: " << e1-s1 << " bp" << std::endl;

				bamSlave.SetRegion( s_frame1.getContigId(), s2, s_frame1.getContigId(), e2 ); //bamSlave.SetRegion( id, s1, id, e2 );
                while( bamSlave.GetNextAlignmentCore(align) )
                {
                    if( !align.IsPaired() ) continue;
                    if( !align.IsMapped() || !align.IsMateMapped() ) continue;

                    if( align.RefID != align.MateRefID ) continue; // pairs must align in the same contig

                    align.BuildCharData(); // fill string fields

                    // se la molteplicità non è stata definita, assumo che sia pari ad 1
                    if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
                    if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';

                    if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

					if( s_frame1.getBegin() <= s_frame2.getBegin() && s2 <= align.Position && align.Position <= e2 && align.MatePosition < s2 ) s_pairs++;
					if( s_frame1.getBegin() > s_frame2.getBegin() && s2 <= align.Position && align.Position <= e2 && align.MatePosition > e2 ) s_pairs++;
                    //if( sr <= align.Position && align.GetEndPosition()-1 <= er && align.MatePosition >= sr && align.MatePosition <= er ) s_pairs++;
					//if( s1 <= align.Position && align.GetEndPosition()-1 <= e1 && align.MatePosition >= s2 && align.MatePosition <= e2 ) s_pairs++;
                    //if( s2 <= align.Position && align.GetEndPosition()-1 <= e2 && align.MatePosition >= s1 && align.MatePosition <= e1 ) s_pairs++;
                }

                if( options.slaveMpBamFiles.size() > 0 && options.slaveMpBamFiles.at( shared.getSlaveFrame().getAssemblyId() ) != "" )
                {
                    bamSlave.Open( options.slaveMpBamFiles.at( shared.getSlaveFrame().getAssemblyId() ) );
                    bamSlave.OpenIndex( options.slaveMpBamFiles.at( shared.getSlaveFrame().getAssemblyId() ) + ".bai");

					bamSlave.SetRegion( s_frame1.getContigId(), s2, s_frame1.getContigId(), e2 );
					while( bamSlave.GetNextAlignmentCore(align) )
					{
						if( !align.IsPaired() ) continue;
						if( !align.IsMapped() || !align.IsMateMapped() ) continue;

						if( align.RefID != align.MateRefID ) continue; // pairs must align in the same contig

						align.BuildCharData(); // fill string fields

						// se la molteplicità non è stata definita, assumo che sia pari ad 1
						if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
						if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';

						if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

						if( s_frame1.getBegin() <= s_frame2.getBegin() && s2 <= align.Position && align.Position <= e2 && align.MatePosition < s2 ) s_mates++;
						if( s_frame1.getBegin() > s_frame2.getBegin() && s2 <= align.Position && align.Position <= e2 && align.MatePosition > e2 ) s_mates++;
						//if( sr <= align.Position && align.GetEndPosition()-1 <= er && align.MatePosition >= sr && align.MatePosition <= er ) s_mates++;
						//if( s1 <= align.Position && align.GetEndPosition()-1 <= e1 && align.MatePosition >= s2 && align.MatePosition <= e2 ) s_mates++;
						//if( s2 <= align.Position && align.GetEndPosition()-1 <= e2 && align.MatePosition >= s1 && align.MatePosition <= e1 ) s_mates++;
					}
				} // end of slave mates number computation

                std::cerr << " slave interval: " << e2-s2 << " bp" << std::endl;

                std::cerr << " master:\tpairs: " << m_pairs << "\tmates: "<< m_mates << std::endl;
				std::cerr << "  slave:\tpairs: " << s_pairs << "\tmates: "<< s_mates << std::endl;
				std::cerr << "----------------------------------------------------" << std::endl;

                bamMaster.Close();
                bamSlave.Close();

//                std::cout << z << std::endl;
//
//                std::set< int32_t > contigs;
//                for(VertexIterator v=vbegin; v!=vend; v++) contigs.insert( ag.getBlock(*v).getMasterFrame().getContigId() );
//                for(std::set<int32_t>::iterator c = contigs.begin(); c != contigs.end(); c++) std::cout << *c << " ";
//                std::cout << std::endl; contigs.clear();
//                for(VertexIterator v=vbegin; v!=vend; v++) contigs.insert( ag.getBlock(*v).getSlaveFrame().getContigId() );
//                for(std::set<int32_t>::iterator c = contigs.begin(); c != contigs.end(); c++) std::cout << *c << " ";
//                std::cout << "\n" << std::endl;
            }*/

            /*ag.removeForks();

            // DEBUG - Assembly graph with forks removed
            std::stringstream ff2;
            ff2 << "./gam_graphs/AssemblyGraph_" << z << "_forks_removed.dot";
            std::ofstream ss2( ff2.str().c_str() );
            ag.writeGraphviz(ss2);
            ss2.close();
            // END of DEBUG */

            if( num_forks == 0 )
            {
				ag_linears++;

                std::stringstream ff2;
                ff2 << "./gam_graphs/AssemblyGraph_" << z << "_linear.dot";
                std::ofstream ss2( ff2.str().c_str() );
                ag.writeGraphviz(ss2);
                ss2.close();

                if( boost::num_vertices(ag) > 0 )
                {
                    //std::vector<Block> filtered = ag.getBlocksVector();
                    //blocksList.push_back(filtered);
                    blocksList.push_back(ag.getBlocksVector());
                }
            }
        }
        catch( boost::not_a_dag ) // if the graph is not a DAG, remove cycles.
        {
			ag_cycles++;

            /* DEBUG - Assembly graph with cycles */
            std::stringstream ff1;
            ff1 << "./gam_graphs/AssemblyGraph_" << z << "_cyclic.dot";
            std::ofstream ss1( ff1.str().c_str() );
            ag.writeGraphviz(ss1);
            ss1.close();
            /* END DEBUG*/

            //ag.removeCycles();
            //cyclesRemoved = true;
        }

        z++; // increase assembly graph counter

        /*if( boost::num_vertices(ag) > 0 )
        {
            std::vector<Block> filtered = ag.getBlocksVector();
            blocksList.push_back(filtered);
        }*/
    }

    _g_statsFile << "[graphs stats]\n"
		<< "Linears = " << ag_linears << "\n"
		<< "Forks = " << ag_forks << "\n"
		<< "Cyclics = " << ag_cycles << "\n"
		<< std::endl;

    return blocksList;
}


std::vector< std::vector< Block > >
partitionBlocksByPairedContigs( const std::vector< Block > &blocks )
{
    typedef PairedContigGraph<> PCGraph;
    typedef boost::graph_traits< PCGraph >::vertex_descriptor Vertex;

    PCGraph pcg(blocks);

    // compute connected components
    std::vector< UIntType > component( boost::num_vertices(pcg) );
    int num = boost::connected_components( pcg, &component[0] );

    // each connected component consists of blocks of contigs which may be extended
    // through weaving.
    std::vector< std::vector< Block > > pairedContigs(num);
    for( UIntType i = 0; i != blocks.size(); i++ )
    {
        Vertex v = pcg.getMasterVertex( blocks[i] );
        pairedContigs.at( component.at(v) ).push_back( blocks[i] );
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
