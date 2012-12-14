/*
 * File:   BuildPctgFunctions.code.hpp
 * Author: vice
 *
 * Created on 4 giugno 2011, 23.29
 */

#include "pctg/BuildPctgFunctions.hpp"
#include "pctg/MergeInCutTailFailed.hpp"

#include "pctg/MergeDescriptor.hpp"

extern MultiBamReader masterBam;
extern MultiBamReader masterMpBam;
extern MultiBamReader slaveBam;
extern MultiBamReader slaveMpBam;

std::list< PairedContig >& buildPctg(
		CompactAssemblyGraph& graph,
		PctgBuilder& builder,
		std::list< PairedContig > &pctgList )
{
	typedef CompactAssemblyGraph::Vertex Vertex;
	typedef CompactAssemblyGraph::VertexIterator VertexIterator;
	typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;
	typedef CompactAssemblyGraph::InEdgeIterator InEdgeIterator;
	
	VertexIterator vbegin,vend;
	MergeBlockLists mergeLists;
	
	// handle forks
	std::vector< MergeBlock > mbv( boost::num_vertices(graph) );
	if( not builder.solveForks(graph,mbv) ) return pctgList;

	// get roots
	std::list< Vertex > roots;
	boost::tie(vbegin,vend) = boost::vertices(graph);
	for( VertexIterator v = vbegin; v != vend; v++ ) if( boost::in_degree(*v,graph) == 0 ) roots.push_back(*v);
	
	// get merge paths (list of merge-block lists)
	while( roots.size() != 0 )
	{
		Vertex rv = roots.back();
		roots.pop_back();
		
		mergeLists.push_front( std::list<MergeBlock>() );
		builder.getMergePaths( graph, rv, mbv, mergeLists );
	}
	
	// DEBUG: print merge list
	/*for( MergeBlockLists::iterator it = mergeLists.begin(); it != mergeLists.end(); it++ )
	{
		for( std::list<MergeBlock>::iterator mb = it->begin(); mb != it->end(); mb++ )
			std::cerr << " (" << mb->m_id << "," << mb->s_id << ")";
		
		std::cerr << std::endl;
	}*/

	for( MergeBlockLists::iterator it = mergeLists.begin(); it != mergeLists.end(); it++ )
		for( std::list<MergeBlock>::iterator mb = it->begin(); mb != it->end(); mb++ )
			builder.alignMergeBlock(graph,*mb);

	builder.splitMergeBlocksByAlign(mergeLists);
	builder.splitMergeBlocksByDirection(mergeLists);
	builder.sortMergeBlocksByDirection(mergeLists);
	builder.splitMergeBlocksByInclusions(mergeLists);

	builder.buildPctgs( pctgList, mergeLists );

	return pctgList;
}


std::list< PairedContig >& buildPctg(
		ThreadedBuildPctg *tbp,
        CompactAssemblyGraph &ag,
        const RefSequence &masterRef,
        const RefSequence &slaveRef,
        std::list< PairedContig > &pctgList )
{
	PctgBuilder builder( tbp, &masterRef, &slaveRef );
	pctgList = buildPctg( ag, builder, pctgList );
	
	return pctgList;
}


void generateSingleCtgPctgs(
        std::list<PairedContig> &pctgList,
        const std::list<int32_t> &ctgIds,
		const RefSequence &masterRef,
        uint64_t &pctgId )
{
	for( std::list<int32_t>::const_iterator i = ctgIds.begin(); i != ctgIds.end(); i++ )
    {
        PctgBuilder builder( NULL, &masterRef, NULL );
        PairedContig pctg = builder.initByContig( pctgId, *i );

        if( pctg.size() > 0 )
        {
            //pctgPool->set( pctg.name(), pctg );
            pctgList.push_back( pctg );
            pctgId++;
        }
    }
}
