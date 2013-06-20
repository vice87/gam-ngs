#include <vector>

#include "strand_fixer/RelativeStrand.hpp"
//#include "boost/graph/graphviz.hpp"

RelativeStrandEvidencesGraph::RelativeStrandEvidencesGraph(const std::list<Block>& blocks) :
        PairedContigGraph<VertexPropType,RelativeStrandEvidences>( blocks )
{
    this->addEdgeWeights(blocks);
}


void
RelativeStrandEvidencesGraph::addEdgeWeights(const std::list<Block>& blocks)
{
    RelativeStrandEvidences evidences;

	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
        Edge e = boost::add_edge(this->getMasterVertex(*b), this->getSlaveVertex(*b), *this).first;
        put( boost::edge_weight_t(), *this, e, evidences );
    }

    for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
        Edge e = boost::edge(this->getMasterVertex(*b), this->getSlaveVertex(*b), *this).first;
        evidences = get( boost::edge_weight_t(), *this, e );

        if( b->getMasterFrame().getStrand() == b->getSlaveFrame().getStrand() )
        {
            evidences.addPositiveEvidences( b->getReadsNumber() );
        }
        else
        {
            evidences.addNegativeEvidences( b->getReadsNumber() );
        }

        put( boost::edge_weight_t(), *this, e, evidences );
    }
}


void
RelativeStrandEvidencesGraph::initVerticesColor()
{
    VertexIterator begin, end;
    boost::tie(begin,end) = boost::vertices(*this);

    for( VertexIterator v = begin; v != end; v++ )
    {
        boost::put( boost::vertex_color_t(), *this, *v, WHITE );
    }
}


void
RelativeStrandEvidencesGraph::extendPathFrom(
        PathsToProbListMap& pathToProb,
        const Vertex& node,
        const RealType& posStrandPathProb,
        UIntType minPathEvidences)
{
    // change node's color to grey
    boost::put(boost::vertex_color_t(), *this, node, GREY);

    AdjacencyIterator begin,end;
    boost::tie(begin,end) = boost::adjacent_vertices(node,*this);

    for( AdjacencyIterator v = begin; v != end; v++ )
    {
        ColorType color = boost::get(boost::vertex_color_t(), *this, *v);

        // if the node hasn't been visited yet
        if( color == WHITE )
        {
            // recover the edge (node,*v) strand evidences
            Edge e = boost::edge(node,*v,*this).first;
            RelativeStrandEvidences edgeStrandProb = boost::get(boost::edge_weight_t(), *this, e);

            // compute the probability of the new path
            RealType newPathProb = this->composePosStrandProb( posStrandPathProb, edgeStrandProb );

            minPathEvidences = std::min( minPathEvidences, edgeStrandProb.getEvidences() );

            // memorize new acyclic path probability
            if( pathToProb[*v].size() < MAX_PTP_LIST_SIZE )
            {
                pathToProb[*v].push_front( EvidenceProbPair(minPathEvidences,newPathProb) );

                // extend path from *v
                this->extendPathFrom( pathToProb, *v, newPathProb, minPathEvidences );
            }
        }
    }

    // we want to compute all the acyclic path, hence
    // all the edges after a visit have to be re-white
    boost::put(boost::vertex_color_t(), *this, node, WHITE);
}


RealType
RelativeStrandEvidencesGraph::composePosStrandProb(
        const RealType& pathPosStrandProb,
        const RelativeStrandEvidences& edgeEvidences)
{
    return pathPosStrandProb * edgeEvidences.getPositiveStrandProb() +
           (1-pathPosStrandProb) * edgeEvidences.getNegativeStrandProb();
}


void
RelativeStrandEvidencesGraph::computePathToProbFrom(PathsToProbListMap& pathToProb, const Vertex& node)
{
    VertexIterator vbegin, vend;
    boost::tie(vbegin,vend) = boost::vertices(*this);

    for(VertexIterator v = vbegin; v != vend; v++ )
    {
        pathToProb[node] = EvidenceProbPairList();
    }

    UIntType evidences(1);
    RealType prob(1);

    pathToProb[node].push_front( EvidenceProbPair(evidences,prob) );

    RealType posStrandPathProb(1);

    // change node's color to grey
    boost::put( boost::vertex_color_t(), *this, node, GREY );

    AdjacencyIterator begin, end;
    boost::tie(begin,end) = boost::adjacent_vertices(node,*this);

    for( AdjacencyIterator v = begin; v != end; v++ )
    {
        ColorType color = boost::get( boost::vertex_color_t(), *this, *v );

        // if the node hasn't been visited yet
        if( color == WHITE )
        {
            // get the edge (node,*v) strand evidences
            Edge e = boost::edge( node, *v, *this ).first;
            RelativeStrandEvidences edgeStrandProb = boost::get( boost::edge_weight_t(), *this, e );

            // compute the probability of the new path
            RealType newPathProb = this->composePosStrandProb( posStrandPathProb, edgeStrandProb );

            // memorize new acyclic path probability
            pathToProb[*v].push_front( EvidenceProbPair(edgeStrandProb.getEvidences(),newPathProb) );

            // extend path from *v
            this->extendPathFrom( pathToProb, *v, newPathProb, edgeStrandProb.getEvidences() );
        }
    }

    // we want to compute all the acyclic path, hence
    // all the edges after a visit have to be re-white
    boost::put( boost::vertex_color_t(), *this, node, WHITE );
}


RelativeStrandEvidencesGraph::StrandProbMapPair
RelativeStrandEvidencesGraph::computeRelativeStrandsWithRespectTo(const Vertex &node)
{
    StrandProbMapPair output;

    if( boost::num_vertices(*this) > 1 )
    {
        this->initVerticesColor();
        PathsToProbListMap pathToProb;
        this->computePathToProbFrom( pathToProb, node );

        PathsToProbListMap::iterator map;
        for( map = pathToProb.begin(); map != pathToProb.end(); map++ )
        {
            Vertex node = map->first;
            RealType vertexProb = 0;
            IntType pathMinEvidences = 0;

            EvidenceProbPairList::iterator prob;
            for( prob = (map->second).begin(); prob != (map->second).end(); prob++ )
            {
                vertexProb += (prob->second * ((RealType)(prob->first)));
                pathMinEvidences += (prob->first);
            }

            vertexProb = vertexProb / ((RealType) pathMinEvidences);

            if( this->isMasterNode(node) )
                output.first[ this->_vertexToCtg[node] ] = vertexProb;
            else
                output.second[ this->_vertexToCtg[node] ] = vertexProb;
        }
    }
    else
    {
        if( this->isMasterNode(node) )
            output.first[ this->_vertexToCtg[node] ] = RealType(1);
        else
            output.second[ this->_vertexToCtg[node] ] = RealType(1);
    }

    return output;
}


RelativeStrandEvidencesGraph::StrandProbMapPair
RelativeStrandEvidencesGraph::computeRelativeStrands()
{
    VertexIterator begin, end;
    boost::tie(begin,end) = boost::vertices(*this);

    if( begin != end ) return this->computeRelativeStrandsWithRespectTo(*begin);

    return StrandProbMapPair();
}


std::pair< std::map<int32_t,StrandProbability>, std::map<int32_t,StrandProbability> >
computeRelativeStrandMap(const std::list<Block>& blocks)
{
    // build strand graph
    RelativeStrandEvidencesGraph rseg(blocks);

    return rseg.computeRelativeStrands();
}
