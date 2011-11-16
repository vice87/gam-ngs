
#include <vector> 
#include <sstream>
#include <fstream>

#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <iostream>

#include "PartitionFunctions.hpp"
#include "graphs/PairedGraph.code.hpp"
#include "graphs/AssemblyGraph.hpp"

std::list< std::vector<Block> >
partitionBlocks( const std::vector<Block> &blocks )
{
    std::list< std::vector<Block> > blocksList;
    
    // group blocks between contigs which may be merged together through weaving
    std::vector< std::vector<Block> > pairedContigsBlocks =
            partitionBlocksByPairedContigs( blocks );
    
    int z = 1; // counter for debug.
    
    // for each partition of blocks
    std::vector< std::vector<Block> >::iterator pcb;
    for( pcb = pairedContigsBlocks.begin(); pcb != pairedContigsBlocks.end(); pcb++ )
    {
        //if(z==55) std::cout << "-----------------------------------------------" << std::endl;
        
        // create an assembly graph
        AssemblyGraph ag( *pcb );
        
        typedef std::vector< size_t > container;
        container c;
        bool done = false;
        bool cyclesRemoved = false;
        
        do
        {
            try
            {
                boost::topological_sort( ag, std::back_inserter(c) );
                done = true;
                
                /* DEBUG - Assembly graph with "fork" nodes */
                typedef boost::graph_traits<AssemblyGraph>::vertex_iterator VertexIterator;
                VertexIterator vbegin,vend;
                boost::tie(vbegin,vend) = boost::vertices(ag);
                
                for (VertexIterator v=vbegin; v!=vend; v++) 
                {
                    int in = boost::in_degree(*v,ag);
                    int out = boost::out_degree(*v,ag);
                    
                    if( in > 1 || out > 1)
                    {
                        std::stringstream ff1;
                        ff1 << "./tmp/Partition_" << z << "_forks.dot";
                        std::ofstream ss1( ff1.str().c_str() );
                        ag.writeGraphviz(ss1);
                        ss1.close();
                        
                        break;
                    }
                }
                /* END DEBUG */
                
                ag.removeForks();
                
                /* DEBUG: OUTPUT ASSEMBLY GRAPH WITHOUT FORKS */
                std::stringstream ff1;
                ff1 << "./tmp/Partition_" << z << "_forks_removed.dot";
                std::ofstream ss1( ff1.str().c_str() );
                ag.writeGraphviz(ss1);
                ss1.close();
                /* END DEBUG*/
            }
            catch( boost::not_a_dag ) // if the graph is not a DAG, remove cycles.
            {
                /* DEBUG - Assembly graph with cycles */
                std::stringstream ff1;
                ff1 << "./tmp/Partition_" << z << "_cycles.dot";
                std::ofstream ss1( ff1.str().c_str() );
                ag.writeGraphviz(ss1);
                ss1.close();
                /* END DEBUG*/
                
                ag.removeCycles();
                cyclesRemoved = true;
            }
        } while(!done);
        
        z++;
        
        if( boost::num_vertices(ag) > 0 )
        {
            std::vector<Block> filtered = ag.getBlocksVector();
            blocksList.push_back(filtered);
        }
    }
    
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
    
    return pairedContigs;
}
