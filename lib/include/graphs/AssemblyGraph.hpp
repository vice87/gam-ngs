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
 * \file AssemblyGraph.hpp
 * \brief Definition of AssemblyGraph class.
 * \details This file contains the definition of the class representing the graph of assemblies.
 */

#ifndef ASSEMBLYGRAPH_HPP
#define	ASSEMBLYGRAPH_HPP

#include <vector>
#include <iostream>
#include <iomanip>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
//#include <boost/graph/graphviz.hpp>

#include "assembly/Block.hpp"
#include "OrderingFunctions.hpp"
#include "strand_fixer/RelativeStrand.hpp"
#include "strand_fixer/StrandProbability.hpp"

namespace boost
{
    enum edge_kind_t { edge_kind };
    BOOST_INSTALL_PROPERTY(edge, kind);
}

typedef enum { MASTER_EDGE, SLAVE_EDGE, BOTH_EDGE } __attribute__((packed)) EdgeKindType;

struct EdgeProperty
{
	EdgeKindType kind;
	double weight;
	int32_t rnum;
	bool min_cov;
};

//! Class implementing the graph of assemblies
/*!
 * The graph is constructed from a vector of blocks. For each block, a relative
 * node is created. Edges connect nodes according to the order of the relative
 * blocks in the master and slave assemblies.
 */
class AssemblyGraph : public boost::adjacency_list< boost::setS, boost::vecS, boost::bidirectionalS,
boost::no_property, boost::property<boost::edge_kind_t,EdgeProperty> >
{
public:
    //typedef boost::adjacency_list< boost::setS, boost::vecS, boost::bidirectionalS,
	//	boost::no_property, boost::property<boost::edge_kind_t,EdgeProperty> > Graph;
    typedef boost::graph_traits<AssemblyGraph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<AssemblyGraph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<AssemblyGraph>::adjacency_iterator AdjacencyIterator;
    typedef boost::graph_traits<AssemblyGraph>::edge_descriptor Edge;
    typedef boost::graph_traits<AssemblyGraph>::edge_iterator EdgeIterator;
    typedef std::map< int32_t, StrandProbability > StrandProbMap;

private:
    uint64_t _agId;						//!< assemblies' graph id
	std::vector<Block> _blockVector;	//!< vector of the blocks used as nodes in the graph

    //! Initialize the graph from a vector of blocks.
    /*!
     * Creates the nodes of the blocks and connects them with a directed edge,
     * according to their orders in master and slave assemblies.
     *
     * \param blocks a vector of block
     */
    void initGraph( const std::list<Block> &blocks );

    //! Adds blocks' vertices to the graph.
    UIntType addVertices();

    //! Connects a block (node) to his successive blocks (nodes) in the master assembly
    /*!
     * \param vertex index of a vertex
     * \param strandMap StranProbMap object
     * \param indexes of the ordered blocks
     * \param back indexes of the ordered blocks
     */
    void addMasterEdges( const UIntType &vertex,
                         const StrandProbMap &strandMap,
                         const std::vector<UIntType> &index,
                         const std::vector<UIntType> &backIndex );

    //! Connects a block (node) to his successive blocks (nodes) in the master assembly
    /*!
     * \param vertex index of a vertex
     * \param strandMap StranProbMap object
     * \param indexes of the ordered blocks
     * \param back indexes of the ordered blocks
     */
    void addSlaveEdges( const UIntType &vertex,
                         const StrandProbMap &strandMap,
                         const std::vector<UIntType> &index,
                         const std::vector<UIntType> &backIndex );

    //! Connects two vertices, whose blocks are successive in the master assebly.
    /*!
     * The edge is directed from \c s to \t
     *
     * \param s first vertex
     * \param t second vertex
     */
    bool addMasterSingleEdge( const UIntType& s, const UIntType& t );

    //! Connects two vertices, whose blocks are successive in the slave assebly.
    /*!
     * The edge is directed from \c s to \t
     *
     * \param s first vertex
     * \param t second vertex
     */
    bool addSlaveSingleEdge( const UIntType& s, const UIntType& t );

    void bubbleDFS( Vertex v, std::vector<char> &colors, bool &found );

public:

    //! A constructor with no arguments.
    /*!
     * Creates an empty graph.
     */
    AssemblyGraph( uint64_t id = 0 );

    //! A constructor.
    /*!
     * Creates a graph of assemblies, given a list of blocks.
     * \param blocks a list of blocks.
     */
	AssemblyGraph( const std::list< Block > &blocks, uint64_t id = 0 );

	inline uint64_t getId() const { return this->_agId; }

    //! Gets the block vector of the graph.
    /*!
     * \return a reference to the block vector
     */
    const std::vector<Block>& getBlocksVector() const;

    //! Gets a block of the graph.
    /*!
     * \param pos index of the block in the vector \c _blockVector
     * \return the block at index \c pos in the block's vector
     */
    const Block& getBlock( const UIntType &pos ) const;

    //! Remove cycles from the graph.
    /*!
     * \return a vector representing the strongly connected components.
     */
    void removeCycles();

    //std::list<Block> removeCyclesFromSCC( std::list<Block> &sccBlocks );

    //! Remove nodes with in/out degree greater than 1, to perform only safe merges.
    void removeForks();

    //! Assign operator of the AssemblyGraph class.
    const AssemblyGraph& operator=( const AssemblyGraph &orig );

    std::ostream& writeGraphviz(std::ostream& os);

    void reverseEdges();

    static void agTopologicalSort( const AssemblyGraph &g, std::list<Vertex> &tsList );
    static void agTopologicalSort( const AssemblyGraph &g, Vertex v, std::vector<char> &colors, std::list<Vertex> &tsList );

	bool hasForks();
	bool hasBubbles();
};

#endif	/* ASSEMBLYGRAPH_HPP */

