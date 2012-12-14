/*!
 * \file CompactAssemblyGraph.hpp
 * \brief Definition of CompactAssemblyGraph class.
 * \details This file contains the definition of the class representing the graph of assemblies collapsing blocks with the same master and slave.
 */

#ifndef COMPACTASSEMBLYGRAPH_HPP
#define	COMPACTASSEMBLYGRAPH_HPP

#include <vector>
#include <list>

#include <boost/graph/adjacency_list.hpp>

#include "assembly/Block.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "strand_fixer/StrandProbability.hpp"

//! Class implementing the graph of assemblies
/*!
 * The graph is constructed from a vector of blocks. For each block, a relative
 * node is created. Edges connect nodes according to the order of the relative
 * blocks in the master and slave assemblies.
 */
class CompactAssemblyGraph : public boost::adjacency_list< boost::setS, boost::vecS,
	boost::bidirectionalS, boost::no_property, boost::property<boost::edge_kind_t,EdgeProperty> >
{
public:
    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::bidirectionalS,
		boost::no_property, boost::property<boost::edge_kind_t,EdgeProperty> > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
	typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
	typedef boost::graph_traits<Graph>::in_edge_iterator InEdgeIterator;

private:
	uint64_t _cgId;
    uint64_t _num_vertices;
	std::vector< std::list<Block> > _blockVector; //!< associate at each vertex a list of blocks (with the same master/slave)

    //! Initialize the graph from a vector of blocks.
    /*!
     * Creates the nodes of the blocks and connects them with a directed edge,
     * according to their orders in master and slave assemblies.
     *
     * \param blocks a vector of block
     */
    void initGraph( const AssemblyGraph &ag );

	void initGraphDFS( const AssemblyGraph &ag, Vertex v, std::vector<char> &colors, std::vector<Vertex> &ag2cg, Vertex u );
	
	bool bubbleDFS( Vertex v, std::vector<char> &colors, bool &found );
	
	std::pair<double,int32_t>
	getRegionScore( MultiBamReader &peBamReader, MultiBamReader &mpBamReader, EdgeKindType kind, 
						   std::list<Block>& b1, std::list<Block>& b2 );
	
	std::pair< std::vector<double>, std::vector<int32_t> >
	getLibRegionScore2( MultiBamReader &bamReader, EdgeKindType kind, std::list<Block>& b1, std::list<Block>& b2 );

public:

    //! A constructor.
    /*!
     * Creates a graph of assemblies, given a list of blocks.
     * \param blocks a list of blocks.
     */
	CompactAssemblyGraph( const AssemblyGraph &ag );
	
	inline uint64_t getId() const { return this->_cgId; }

    //! Gets the blocks vector of the graph.
    /*!
     * \return a reference to the blocks vector
     */
    const std::vector< std::list<Block> >& getBlocksVector() const;

    //! Gets a block of the graph.
    /*!
     * \param pos index of the block in the vector \c _blockVector
     * \return the block at index \c pos in the block's vector
     */
    const std::list<Block>& getBlocks( const Vertex &pos ) const;

    //! Assign operator of the AssemblyGraph class.
    const CompactAssemblyGraph& operator=( const CompactAssemblyGraph &orig );
	
	bool hasBubbles();

    std::ostream& writeGraphviz(std::ostream& os);
	
	void computeEdgeWeights( MultiBamReader &masterBamReader, MultiBamReader &masterMpBamReader, 
							 MultiBamReader &slaveBamReader, MultiBamReader &slaveMpBamReader );
};

#endif	/* COMPACTASSEMBLYGRAPH_HPP */

