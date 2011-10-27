/*! 
 * \file AssemblyGraph.hpp
 * \brief Definition of AssemblyGraph class.
 * \details This file contains the definition of the class representing the graph of assemblies.
 */

#ifndef ASSEMBLYGRAPH_HPP
#define	ASSEMBLYGRAPH_HPP

#include <boost/graph/adjacency_list.hpp>

#include "assembly/Block.hpp"
#include "strand_fixer/StrandProbability.hpp"

namespace boost 
{
    enum edge_kind_t { edge_kind };
    BOOST_INSTALL_PROPERTY(edge, kind);
}

typedef enum { MASTER_EDGE, SLAVE_EDGE, BOTH_EDGE } EdgeKindType;

//! Class implementing the graph of assemblies
/*!
 * The graph is constructed from a vector of blocks. For each block, a relative
 * node is created. Edges connect nodes according to the order of the relative
 * blocks in the master and slave assemblies.
 */
class AssemblyGraph : public boost::adjacency_list< boost::setS, boost::vecS, 
        boost::bidirectionalS, boost::no_property, boost::property<boost::edge_kind_t,EdgeKindType> >
{
public:
    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::bidirectionalS, 
            boost::no_property, boost::property<boost::edge_kind_t,EdgeKindType> > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef std::map< IdType, StrandProbability > StrandProbMap;
    
private:
    std::vector<Block> _blockVector; //!< vector of the blocks used as nodes in the graph
    
    //! Initialize the graph from a vector of blocks.
    /*!
     * Creates the nodes of the blocks and connects them with a directed edge,
     * according to their orders in master and slave assemblies.
     * 
     * \param blocks a vector of block
     */
    void initGraph( const std::vector<Block> &blocks );
    
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


public:
    
    //! A constructor with no arguments.
    /*!
     * Creates an empty graph.
     */
    AssemblyGraph();
    
    //! A constructor.
    /*!
     * Creates a graph of assemblies, given a vector of blocks.
     * \param blocks a vector of blocks.
     */
    AssemblyGraph( const std::vector< Block > &blocks );
    
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
    
    std::list<Block> removeCyclesFromSCC( std::list<Block> &sccBlocks );
    
    //! Assign operator of the AssemblyGraph class.
    const AssemblyGraph& operator=( const AssemblyGraph &orig );
    
    std::ostream& writeGraphviz(std::ostream& os);

};

#endif	/* ASSEMBLYGRAPH_HPP */

