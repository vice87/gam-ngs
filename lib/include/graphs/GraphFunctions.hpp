/*! 
 * \file GraphFunctions.hpp
 * \brief Definition of utility functions for graphs.
 */

#ifndef _GRAPH_FUNCTIONS_
#define _GRAPH_FUNCTIONS_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

//! Returns a graph with reversed edge.
/*!
 * \param g a graph.
 * \return the graph \c with reversed edge
 */
template <class GRAPH>
inline GRAPH
get_reversed_edges_graph(const GRAPH& g)
{
  GRAPH out(g);
  
  typename boost::graph_traits<GRAPH>::edge_iterator e_i, e_end;
  typename boost::graph_traits<GRAPH>::edge_descriptor e;
  bool found;

  boost::tie(e_i,e_end)=boost::edges(g);

  while (e_i!=e_end) {
    boost::add_edge( boost::target(*e_i,g), boost::source(*e_i,g), out );
    boost::tie(e,found) = boost::edge( boost::source(*e_i,g), boost::target(*e_i,g), out );

    if (found) boost::remove_edge(e,out);

    e_i++;
  }

  return out;
}

//! Returns whether a graph contains vertices with out degree greater than 1.
/*!
 * \param g a graph.
 * \return \c true if the maximum out degree of all vertices is 1.
 */
template <class GRAPH>
inline bool
graph_forks(const GRAPH& g)
{
  typename boost::graph_traits<GRAPH>::vertex_iterator v_i, v_end;

  boost::tie(v_i,v_end)=boost::vertices(g);

  while (v_i!=v_end) 
  {
    if (boost::out_degree(*v_i, g) > 1) return true;

    v_i++;
  }

  return false;
}

//! Returns whether a graph contains vertices with out/in degree greater than 1.
/*!
 * \param g a graph.
 * \return \c true if the maximum in/out degree of all vertices is 1.
 */
template <class GRAPH>
inline bool
is_linear(const GRAPH& g)
{
  if( graph_forks(g) || graph_forks(get_reversed_edges_graph(g)) ) 
      return false;

  return true;
}

#endif // _GRAPH_FUNCTIONS_
