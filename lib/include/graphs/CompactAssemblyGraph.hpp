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

/*!
 * \file CompactAssemblyGraph.hpp
 * \brief Definition of CompactAssemblyGraph class.
 * \details This file contains the definition of the class representing the graph of assemblies collapsing blocks with the same master and slave.
 */

#ifndef COMPACTASSEMBLYGRAPH_HPP
#define	COMPACTASSEMBLYGRAPH_HPP

#include <vector>
#include <list>
#include <iostream>
#include <iomanip>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
//#include <boost/graph/graphviz.hpp>
#include <boost/dynamic_bitset.hpp>

#include "OrderingFunctions.hpp"
#include "assembly/Block.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "strand_fixer/RelativeStrand.hpp"
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
//    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::bidirectionalS,
//		boost::no_property, boost::property<boost::edge_kind_t,EdgeProperty> > Graph;
    typedef boost::graph_traits<CompactAssemblyGraph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<CompactAssemblyGraph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<CompactAssemblyGraph>::adjacency_iterator AdjacencyIterator;
    typedef boost::graph_traits<CompactAssemblyGraph>::edge_descriptor Edge;
    typedef boost::graph_traits<CompactAssemblyGraph>::edge_iterator EdgeIterator;
    typedef boost::graph_traits<CompactAssemblyGraph>::out_edge_iterator OutEdgeIterator;
    typedef boost::graph_traits<CompactAssemblyGraph>::in_edge_iterator InEdgeIterator;

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
	void initGraph2( const AssemblyGraph &ag );

	void initGraphDFS_NR(
		const AssemblyGraph &ag,
		const AssemblyGraph::Vertex &root,
		boost::dynamic_bitset<> *visited,
		std::vector<Vertex> *ag2cg
	);

    void initGraph( const AssemblyGraph &ag );

    void initGraphDFS(
        const AssemblyGraph &ag,
        const AssemblyGraph::Vertex &v,
		const AssemblyGraph::Vertex &u,
        boost::dynamic_bitset<> *colors,
        std::vector<Vertex> *ag2cg
	);

    void bubbleDFS( Vertex v, std::vector<char> &colors, bool &found );

    void getRegionScore( MultiBamReader &peBamReader, MultiBamReader &mpBamReader, EdgeKindType kind,
            std::list<Block>& b1, std::list<Block>& b2, double &weight, int32_t &rnum, bool &min_cov );

    void getLibRegionScore( MultiBamReader &bamReader, EdgeKindType kind, std::list<Block> &b1, std::list<Block> &b2,
            double &weight, int32_t &rnum, bool &min_cov );

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

