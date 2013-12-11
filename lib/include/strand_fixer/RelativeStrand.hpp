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
 * \file RelativeStrand.hpp
 * \brief Definition of RelativeStrandEvidencesGraph class.
 * \details This file contains the definition of a graph built on blocks that
 * allows to estimate the probability of a contig to be reverse complemented.
 */

#ifndef RELATIVESTRAND_HPP
#define	RELATIVESTRAND_HPP

#include <map>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "graphs/PairedGraph.code.hpp"
#include "strand_fixer/StrandProbability.hpp"
#include "strand_fixer/RelativeStrandEvidences.hpp"

#define MAX_PTP_LIST_SIZE 100

//! Class implementing a relative strand graph.
class RelativeStrandEvidencesGraph : public PairedContigGraph< unsigned short int, RelativeStrandEvidences >
{
    //! Enum representing the color of a node in the graph.
    enum Colors
    {
        WHITE,
        GREY,
        BLACK
    } __attribute__((packed));

    typedef boost::property<boost::edge_weight_t, RelativeStrandEvidences> EdgeStrandProbabilityProperty;
    typedef RelativeStrandEvidences EdgeWeightType;
    typedef unsigned short int VertexPropType;

    typedef boost::adjacency_list< boost::setS, boost::vecS, boost::undirectedS,
                                   boost::property<boost::vertex_color_t, VertexPropType >,
                                   boost::property<boost::edge_weight_t, EdgeWeightType> > Graph;

    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
    typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;

    typedef PairedGraph< VertexPropType, EdgeWeightType > PairedGraphType;

    typedef std::pair< UIntType, RealType > EvidenceProbPair;
    typedef std::list<EvidenceProbPair> EvidenceProbPairList;
    typedef std::map< Vertex, EvidenceProbPairList > PathsToProbListMap;

    typedef std::map< int32_t, StrandProbability > StrandProbMap;
    typedef std::pair<StrandProbMap, StrandProbMap> StrandProbMapPair;

    typedef VertexPropType ColorType;

    typedef std::map< IdType, unsigned int > LabelMapType;
    typedef std::pair< LabelMapType, LabelMapType> FullLabelMapType;
    typedef std::vector<IdType> VertexLabelType;

    void initVerticesColor(); //

    void extendPathFrom
    (
        PathsToProbListMap &pathToProb,
        const Vertex &node,
        const RealType &posStrandPathProb,
        UIntType minPathEvidences
    );

    RealType composePosStrandProb
    (
        const RealType &pathPosStrandProb,
        const RelativeStrandEvidences &edgeEvidences
    );

    void computePathToProbFrom //
    (
        PathsToProbListMap &pathToProb,
        const Vertex &node
    );

    void addEdgeWeights( const std::list<Block> &blocks ); //

public:

    RelativeStrandEvidencesGraph(const std::list<Block> &blocks); //

    StrandProbMapPair computeRelativeStrandsWithRespectTo(const Vertex &node); //

    StrandProbMapPair computeRelativeStrands(); //

}; // class RelativeStrandEvidencesGraph


std::pair< std::map<int32_t,StrandProbability>, std::map<int32_t,StrandProbability> >
computeRelativeStrandMap( const std::list<Block> &blocks );

#endif	/* RELATIVESTRAND_HPP */

