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
    std::list< PairedContig > pctgList; // add
    PairedContig pctg(pctgId);
    
    if( boost::num_vertices(AG) == 0 ) return pctgList;
    
    //typedef std::vector< UIntType > container;
    //container c;
    
    typedef boost::graph_traits<AssemblyGraph>::vertex_descriptor Vertex;
    std::list< Vertex > c;
    std::list< UIntType > notMerged;
    
    try
    {
	AssemblyGraph::agTopologicalSort( AG, c ); //boost::topological_sort(AG,std::back_inserter(c));
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
        //std::vector<UIntType>::reverse_iterator i;
        std::list<Vertex>::reverse_iterator i;
        
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
                        pctg = builder.extendByBlock( pctg, AG.getBlock(notMerged.front()), AG.getBlock(notMerged.back()) );
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
                pctg = builder.extendByBlock( pctg, AG.getBlock(notMerged.front()), AG.getBlock(notMerged.back()) );
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

