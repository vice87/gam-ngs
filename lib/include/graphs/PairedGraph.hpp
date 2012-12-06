/*!
 * \file PairedGraph.hpp
 * \brief Definition of PairedGraph and PairedContigGraph classes.
 * \details These classes allow to construct an undirected graph of contigs of
 * master and slave assemblies which are linked toghether iff there is a block
 * on them.
 */

#ifndef PAIRED_GRAPH_
#define PAIRED_GRAPH_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "types.hpp"
#include "assembly/Block.hpp"

//! PairedGraph class.
/*!
 * Base class for paired graphs.
 */
template <class VERTEX_PROP = UIntType, class EDGE_WEIGHT = UIntType>
class PairedGraph : public boost::adjacency_list< boost::setS, boost::vecS, boost::undirectedS,
        boost::property< boost::vertex_color_t, VERTEX_PROP >,
        boost::property< boost::edge_weight_t, EDGE_WEIGHT > >
{

public:
    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::undirectedS,
                boost::property< boost::vertex_color_t, VERTEX_PROP >,
                boost::property< boost::edge_weight_t, EDGE_WEIGHT > > Graph; //!< Graph type.

    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge; //!< Edge descriptor type.
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex; //!< Vertex descriptor type.

protected:
    Vertex _firstSlaveVertex;                               //!< vertex descriptor of the first slave vertex.
    std::map< int32_t,Vertex > _masterMap; //!< associates a master contig to its vertex.
    std::map< int32_t,Vertex > _slaveMap;  //!< associates a slave contig to its vertex.
    std::vector< int32_t > _vertexToCtg;   //!< associates a vertex to its contig identifier.

public:
    //! A constructor.
    /*!
     * Creates an empty PairedGraph.
     */
    PairedGraph();

    //! A copy constructor.
    /*!
     * Creates a copy of a PairedGraph.
     * \param orig a PairedGraph.
     */
    PairedGraph( const PairedGraph &orig );

    //! Gets the vertex descriptor of the master contig of a block.
    /*!
     * \param block a block.
     * \return vertex descriptor of the master contig of \c block.
     */
    Vertex getMasterVertex( const Block &block ) const;

    //! Gets the vertex descriptor of the slave contig of a block.
    /*!
     * \param block a block.
     * \return vertex descriptor of the slave contig of \c block.
     */
    Vertex getSlaveVertex( const Block &block ) const;

	//! Returns whether a vertex corresponds to a master contig.
    /*!
     * \param node a vertex descriptor.
     * \return \c true if \c node corresponds to a master contig, \c false otherwise.
     */
    bool isMasterNode( const Vertex &node ) const;

    //! Returns whether a vertex corresponds to a slave contig.
    /*!
     * \param node a vertex descriptor.
     * \return \c true if \c node corresponds to a slave contig, \c false otherwise.
     */
    bool isSlaveNode( const Vertex &node ) const;

    //! Assign operator of the PairedContig graph.
    const PairedGraph& operator=( const PairedGraph &orig );
};


//! PairedContigGraph class.
/*!
 * Extends PairedContig class. Connects vertices associated to master and slave contigs of a block.
 */
template <class VERTEX_PROP = UIntType, class EDGE_WEIGHT = UIntType>
class PairedContigGraph : public PairedGraph<VERTEX_PROP,EDGE_WEIGHT>
{

public:
    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::undirectedS,
                boost::property< boost::vertex_color_t, VERTEX_PROP >,
                boost::property< boost::edge_weight_t, EDGE_WEIGHT > > Graph;

    typedef typename PairedGraph<VERTEX_PROP,EDGE_WEIGHT>::Edge Edge; //!< Edge descriptor type.

    typedef typename boost::graph_traits< Graph >::vertex_iterator VertexIterator;
    typedef typename boost::graph_traits< Graph >::edge_iterator EdgeIterator;

private:
    //! Initialises the graph.
    /*!
     * Adds an edge between vertices whose contigs are the master and the slave contig of a block.
     *
     * \param blocks a vector of blocks.
     */
    void initGraph( const std::list<Block>& blocks );

    //! Creates a vertex for each contig.
    /*!
     * \param blocks a vector of blocks.
     */
    void initVertexLabels( const std::list<Block>& blocks );

public:
    //! A constructor.
    /*!
     * Creates an empty PairedContigGraph.
     */
    PairedContigGraph();

    //! A copy constructor.
    /*!
     * Creates a copy of a PairedContigGraph.
     * \param orig a PairedContigGraph object
     */
    PairedContigGraph( const PairedContigGraph &orig );

    //! A constructor.
    /*!
     * Creates and initialises a PairedContigGraph, given a vector of blocks.
     * \param blocks a vector of blocks.
     */
    PairedContigGraph( const std::list<Block> &blocks );

    //! Write the graph in dot format.
    /*!
     * \param os output stream
     */
    std::ostream& writeGraphviz(std::ostream& os);
};

#endif // PAIRED_GRAPH_
