/* 
 * File:   BuildPctgFunctions.code.hpp
 * Author: vice
 *
 * Created on 4 giugno 2011, 23.29
 */

#ifndef BUILDPCTGFUNCTIONS_CODE_HPP
#define	BUILDPCTGFUNCTIONS_CODE_HPP

#include "BuildPctgFunctions.hpp"

std::list< PairedContig > buildPctg(
        IdType& pctgId,
        ThreadedBuildPctg *tbp,
        const AssemblyGraph& AG, 
        const PctgBuilder& builder )
{
    std::list< PairedContig > pctgList; //add
    PairedContig pctg(pctgId);
    
    if( boost::num_vertices(AG) == 0 ) return pctgList;
    
    typedef std::vector< UIntType > container;
    container c;//, c_new;
    std::list< UIntType > notMerged;
    
    try
    {
		boost::topological_sort(AG,std::back_inserter(c));
		
		/* new code: weaving sul longest path */
		/*IntType pathLength[ boost::num_vertices(AG) ];
		IntType prev[ boost::num_vertices(AG) ];
		bool done[ boost::num_vertices(AG) ];
		
		for( int i=0; i < boost::num_vertices(AG); i++ )
		{
			pathLength[i] = 0;
			prev[i] = -1;
			done[i] = false;
		}
		
		typedef boost::graph_traits<AssemblyGraph>::adjacency_iterator AdjacencyIterator;
		AdjacencyIterator begin, end;
		
		for( container::reverse_iterator i = c.rbegin(); i != c.rend(); ++i )
		{
			boost::tie(begin,end) = boost::adjacent_vertices(*i,AG);
			for( AdjacencyIterator v = begin; v != end; v++ )
			{
				if( pathLength[*v] <= pathLength[*i] + AG.getBlock(*v).getReadsNumber() )
				{
					pathLength[*v] = pathLength[*i] + AG.getBlock(*v).getReadsNumber();
					prev[*v] = *i; // update back-pointer of the path
				}
			}
		}
		
		IntType longest = -1; int path_node = -1;
		
		// find longest path end
		for( int i=0; i < boost::num_vertices(AG); i++ )
		{
			if(!done[i] && pathLength[i] > longest)
			{
				longest = pathLength[i];
				path_node = i;
			}
		}
		
		while( longest != -1 )
		{
			int cur_node = path_node;
			
			while(cur_node != -1 && !done[cur_node])
			{
				done[cur_node] = true;
				c_new.push_back(cur_node);
				cur_node = prev[cur_node];
			}
			
			// find longest path end
			longest = -1;
			for( int i=0; i < boost::num_vertices(AG); i++ )
			{
				if(!done[i] && pathLength[i] > longest)
				{
					longest = pathLength[i];
					path_node = i;
				}
			}
		}*/
		/* end new code */
    }
    catch( boost::not_a_dag )
    {
        AssemblyGraph newAG( AG );
        
        std::cerr << "With cycles" << std::endl;
        newAG.removeCycles();
        std::cerr << "Without cycles" << std::endl;
        //std::cerr << newAG.getBlocksVector() << std::endl; // FIX
        
        throw boost::not_a_dag();
    }
    
    UIntType last = *(c.rbegin());
    
    try
    {
        std::vector<UIntType>::reverse_iterator i;
        
        i = c.rbegin();
        while( i != c.rend() )
        {
        //for( i = c.rbegin(); i != c.rend(); i++ )
        //{
            if( !Block::shareContig(AG.getBlock(*i),AG.getBlock(last)) ) // modifica rispetto originale andrebbe shareContig (sia master che slave)
            {
                //std::cout << "not share contig: " << AG.getBlock(*i) << std::endl;
                while( notMerged.size() != 0 )
                {
                    try
                    {
                        //std::cout << "Extending pctg with " << AG.getBlock(notMerged.front()).getMasterFrame().getContigId() << "-" << AG.getBlock(notMerged.front()).getSlaveFrame().getContigId() << std::endl << std::flush;
                        pctg = builder.extendByBlock(pctg,AG.getBlock(notMerged.front()));
                        notMerged.clear();
                    }
                    catch( std::logic_error &e )
                    {
                        //std::cout << "  => Logic error exception1: " << e.what() << std::endl << std::flush;
                        
                        //add
                        pctgList.push_back( pctg );
                        //pctgId = tbp->readPctgNumAndIncrease();
                        pctg = PairedContig();
                        continue;
                        //add
                        
                        notMerged.pop_front();
                    }
                }
            }
            
            //std::cout << "adding to notMerged: " << AG.getBlock(*i).getMasterFrame().getContigId() << "-" << AG.getBlock(*i).getSlaveFrame().getContigId() << std::endl << std::flush;
            
            notMerged.push_front(*i);
            last = *i;
            
            i++;
        }
        
        
        while( notMerged.size() != 0 )
        {
            try
            {
                //std::cout << "Extending (last) pctg with " << AG.getBlock(notMerged.front()).getMasterFrame().getContigId() << "-" << AG.getBlock(notMerged.front()).getSlaveFrame().getContigId() << std::endl << std::flush;
                pctg = builder.extendByBlock(pctg,AG.getBlock(notMerged.front()));
                notMerged.clear();
            }
            catch( std::logic_error &e )
            {
                //std::cout << "  => Logic error exception2: " << e.what() << std::endl;
                notMerged.pop_front();
            }
        }
        
        pctgList.push_back( pctg );
        
    }
    catch(...)
    {
        //return PairedContig(pctgId);
    }
    
    return pctgList;
}


std::list< PairedContig > buildPctg( 
        const AssemblyGraph &ag,
        const HashContigMemPool *masterPool, 
        const HashContigMemPool *slavePool, 
        const BamTools::RefVector *masterRefVector,
        const BamTools::RefVector *slaveRefVector,
        IdType &pctgId,
        ThreadedBuildPctg *tbp )
{
    PctgBuilder builder( masterPool, slavePool, masterRefVector, slaveRefVector );
    
    //std::cout << "Building pctg " << pctgId << std::endl << std::flush;
    
    std::list< PairedContig > pctgList;
    pctgList = buildPctg( pctgId, tbp, ag, builder );
    
    return pctgList;
}


void generateSingleCtgPctgs(
        std::list<PairedContig>& pctgList, 
        const std::list<IdType>& ctgIds, 
        HashContigMemPool* pctgPool, 
        HashContigMemPool* masterPool,
        BamTools::RefVector *masterRefVector,
        IdType& pctgId )
{
    std::list<IdType>::const_iterator i;
    for( i = ctgIds.begin(); i != ctgIds.end(); i++ )
    {
        PctgBuilder builder( masterPool, NULL, masterRefVector, NULL );
        PairedContig pctg = builder.initByContig( pctgId, *i);
        
        if( pctg.size() > 0 )
        {
            //pctgPool->set( pctg.name(), pctg );
            pctgList.push_back( pctg );
            pctgId++;
        }
    }
}

#endif	/* BUILDPCTGFUNCTIONS_CODE_HPP */

