#include <iostream>
#include <exception>

#include "pctg/PctgBuilder.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/smith_waterman.hpp"


PctgBuilder::PctgBuilder(
        const HashContigMemPool *masterPool, 
        const HashContigMemPool *slavePool,
        const BamTools::RefVector *masterRefVector,
        const BamTools::RefVector *slaveRefVector) :
                _masterPool(masterPool), _slavePool(slavePool),
                _masterRefVector(masterRefVector), _slaveRefVector(slaveRefVector),
                _maxAlignment(DEFAULT_MAX_SEARCHED_ALIGNMENT),
                _maxPctgGap(DEFAULT_MAX_GAPS),
                _maxCtgGap(DEFAULT_MAX_GAPS) {}


Contig PctgBuilder::loadMasterContig(const IdType& ctgId) const
{
    std::string refName = this->_masterRefVector->at(ctgId).RefName;
    return (this->_masterPool)->get( refName );
}

Contig PctgBuilder::loadSlaveContig(const IdType& ctgId) const
{
    std::string refName = this->_slaveRefVector->at(ctgId).RefName;
    return (this->_slavePool)->get( refName );
}


PairedContig PctgBuilder::initByContig(const IdType& pctgId, const IdType& ctgId) const
{
    PairedContig pctg(pctgId);
    return this->addFirstContigTo(pctg,ctgId);
}


const PairedContig& PctgBuilder::extendPctgWithCtgFrom(
        PairedContig& orig, 
        const Contig& ctg, 
        ContigInPctgInfo& ctgInfo, 
        const std::pair<UIntType,UIntType>& pos, 
        const std::pair<UIntType,UIntType>& gaps, 
        bool isMasterCtg) const
{
    ctgInfo.setGaps((UIntType)gaps.first - (UIntType)gaps.second);
    
    if(isMasterCtg) orig.getMasterCtgMap()[ctgInfo.getId()] = ctgInfo;
    else orig.getSlaveCtgMap()[ctgInfo.getId()] = ctgInfo;
    
    UIntType pctgPos = pos.first;
    UIntType ctgPos = pos.second;
    
    UIntType lastBases = (orig.size()-pctgPos);
    UIntType newBases = (ctg.size()-ctgPos);
    
    if( lastBases < newBases ) orig.resize( pctgPos + newBases );
    
    for( UIntType i = 0; i < newBases; i++ ) // ORIGINALE => for( UIntType i = ctgPos; i < ctg.size(); i++ )
    {
        // avoid to write an 'N' to the paired contig
        if( i < lastBases && ctg.at(ctgPos + i).base() == N ) continue;
        
        orig.at(pctgPos + i) = ctg.at(ctgPos + i);
        //orig.qual(pctgPos+i) = ctg.qual(i);
    }
    
    return orig;
}


const PairedContig& PctgBuilder::extendPctgWithCtgUpto(
        PairedContig& orig, 
        const Contig& ctg, 
        const ContigInPctgInfo& ctgInfo, 
        const UIntType pctgShift, 
        bool isMasterCtg) const
{
    orig = this->shiftPctgOf(orig,pctgShift);
    
    if(isMasterCtg) orig.getMasterCtgMap()[ctgInfo.getId()] = ctgInfo;
    else orig.getSlaveCtgMap()[ctgInfo.getId()] = ctgInfo;
       
    for( UIntType i = 0; i < pctgShift; i++ )
    {
        orig.at(i) = ctg.at(i);
        //orig.qual(i) = ctg.qual(i);
    }
    
    return orig;
}


PairedContig PctgBuilder::shiftPctgOf(const PairedContig& orig, const UIntType shiftSize) const
{
    Contig contigOut(orig.size() + shiftSize);
    
    for( UIntType i=0; i < orig.size(); i++ )
    {
        contigOut.at(shiftSize + i) = orig.at(i);
        //contigOut.qual(shiftSize + i) = orig.qual(i);
    }
    
    return PairedContig( shiftOf(orig,shiftSize), contigOut );
}



PairedContig PctgBuilder::addFirstContigTo(const PairedContig& pctg, const IdType &ctgId) const
{
    if( pctg.size() != 0 ) throw std::invalid_argument( "Paired contig is not empty" );
    PairedContig out(pctg.getId());
    Contig ctg = this->loadMasterContig(ctgId);
    ContigInPctgInfo ctgInfo(ctgId,ctg.size(),0);
    
    return this->extendPctgWithCtgFrom(out,ctg,ctgInfo,std::pair<UIntType,UIntType>(0,0),
            std::pair<UIntType,UIntType>(0,0),true);
}


PairedContig PctgBuilder::addFirstBlockTo(PairedContig pctg, const Block& block) const
{
    IdType masterCtgId = block.getMasterFrame().getContigId();
    pctg = addFirstContigTo(pctg,masterCtgId);
    
    return this->extendByBlock(pctg,block);
}


PairedContig PctgBuilder::extendByBlock(const PairedContig& pctg, const Block& block) const
{
    if( pctg.size() == 0 ) return this->addFirstBlockTo(pctg,block);
    
    IdType masterCtgId, slaveCtgId;
    
    masterCtgId = block.getMasterFrame().getContigId();
    slaveCtgId  = block.getSlaveFrame().getContigId();
    
    
    if( !pctg.containsMasterCtg(masterCtgId) )
    {
        if( !pctg.containsSlaveCtg(slaveCtgId) )
        {
            throw std::logic_error("Paired contig cannot be extended by this block.");
        }
        else
        {
            // merge the master contig of the block with the paired contig
            return this->mergeContig(pctg,block,true);
        }
    }
    else
    {
        if( !pctg.containsSlaveCtg(slaveCtgId) )
        {
            // merge the slave contig of the block with the paired contig
            return this->mergeContig(pctg,block,false);
        }
    }
    
    return pctg;
}


PairedContig PctgBuilder::mergeContig(PairedContig pctg, const Block& block, bool mergeMasterCtg) const
{
    IdType firstCtgId,  // identifier of the contig to be merged
           secondCtgId; // identifier of the contig already inside the paired contig
    
    if( mergeMasterCtg )
    {
        firstCtgId = block.getMasterFrame().getContigId();
        secondCtgId = block.getSlaveFrame().getContigId();
    }
    else
    {
        firstCtgId = block.getSlaveFrame().getContigId();
        secondCtgId = block.getMasterFrame().getContigId();
    }
    
    UIntType pctgPos; // position in pctg from which the alignment should be searched
    
    if( pctg.getContigInfo(secondCtgId,!mergeMasterCtg).isReversed() )
    {
        UIntType endInCtg = (mergeMasterCtg) ? block.getSlaveFrame().getEnd() : block.getMasterFrame().getEnd();
        pctgPos = pctg.getBasePosition(secondCtgId, endInCtg, !mergeMasterCtg);
    }
    else
    {
        UIntType beginInCtg = (mergeMasterCtg) ? block.getSlaveFrame().getBegin() : block.getMasterFrame().getBegin();
        pctgPos = pctg.getBasePosition(secondCtgId, beginInCtg, !mergeMasterCtg);
    }
    
    // load contig that should be merged with pctg.
    Contig ctg = (mergeMasterCtg) ? this->loadMasterContig(firstCtgId) : this->loadSlaveContig(firstCtgId);    
    Frame frame = (mergeMasterCtg) ? block.getMasterFrame() : block.getSlaveFrame();
    
    // find best alignment between ctg and pctg.
    BestPctgCtgAlignment bestAlign( this->findBestAlignment(pctg,pctg.getContigInfo(secondCtgId,!mergeMasterCtg),
                                                pctgPos, ctg, frame) );
    
    // if the alignment is not good enough, return the paired contig
    if( bestAlign.getAlignment().homology() < MIN_HOMOLOGY ||
        bestAlign.getAlignment().length() < MIN_ALIGNMENT  ||
        bestAlign.getAlignment().length() < ((double)ctg.size())*MIN_ALIGNMENT_QUOTIENT )
    {
        //std::cout << "alignment sucks!" << std::endl << std::flush;
        return pctg;
    }
    
    // avoid to overlap primary assembly contigs
    if( mergeMasterCtg )
    {
        if( sameAssemblyCtgsOverlapedBy(pctg,ctg,bestAlign.getAlignment().b_position_in_a(), mergeMasterCtg) )
        {
            //std::cout << "No overlap for primary assembly constraint disattended." << std::endl << std::flush;
            return pctg;
            throw ConstraintsDisattended("No overlap for primary assembly constraint disattended.");
        }
    }
    
    // merge ctg inside pctg, using alignment informations.
    return this->mergeCtgInPos(pctg,ctg,firstCtgId,bestAlign,mergeMasterCtg);
}


const PairedContig& PctgBuilder::mergeCtgInPos(
        PairedContig& pctg, 
        const Contig& ctg, 
        const IdType& ctgId, 
        const BestPctgCtgAlignment& bestAlign, 
        bool mergeMaster) const
{
    if( mergeMaster ) return this->mergeMasterCtgInPos(pctg,ctg,ctgId,bestAlign);
    
    return this->mergeSlaveCtgInPos(pctg,ctg,ctgId,bestAlign);
}


const PairedContig& PctgBuilder::mergeMasterCtgInPos(
        PairedContig& pctg, 
        const Contig& ctg, 
        const IdType& ctgId, 
        const BestPctgCtgAlignment& bestAlign) const
{
    UIntType pctgShift = 0;
    
    if( bestAlign.getAlignment().b_position_in_a() < 0 )
        pctgShift = (-bestAlign.getAlignment().b_position_in_a());
    
    std::pair<UIntType,UIntType> pos;
    pos = first_match_pos_in( bestAlign.getAlignment() );
    
    ContigInPctgInfo ctgInfo( ctgId, bestAlign );
    
    if(pctgShift > 0)
    {
        ctgInfo.setPosition(0);
        pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, pctgShift, true );
    }
    
    return this->extendPctgWithCtgFrom( pctg,
                        ctg, ctgInfo, pos, std::pair<UIntType,UIntType>(0,0), true );
    
}


const PairedContig& PctgBuilder::mergeSlaveCtgInPos(
        PairedContig& pctg, 
        const Contig& ctg, 
        const IdType& ctgId, 
        const BestPctgCtgAlignment& bestAlign) const
{
    UIntType pctgPos, ctgPos=0, pctgShift = 0, prevSize = pctg.size();
    
    if( bestAlign.getAlignment().b_position_in_a() < 0 )
        pctgShift = (-bestAlign.getAlignment().b_position_in_a());
    
    std::pair<UIntType,UIntType> pos,gaps;
    pos = last_match_pos_in( bestAlign.getAlignment() ); // alignment.hpp
    gaps = gaps_before_last_match_in( bestAlign.getAlignment() );
    
    pctgPos = pos.first + pctgShift;
    ctgPos = pos.second;

    ContigInPctgInfo ctgInfo( ctgId, bestAlign );
    
    if(pctgShift > 0)
    {
        ctgInfo.setPosition(0);
        pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, pctgShift, false );
    }
    
    pos.first = pctgPos;
    
    if( bestAlign.getAlignment().end_b_in_a() <= (UIntType)prevSize )
        pos.second = ctg.size();
    
    return this->extendPctgWithCtgFrom( pctg, ctg, ctgInfo, pos, gaps, false );
    
}


BestPctgCtgAlignment PctgBuilder::findBestAlignment(
        const PairedContig& pctg, 
        const ContigInPctgInfo& pctgInfo, 
        const UIntType pctgPos, Contig& ctg, 
        const Frame ctgFrame) const
{
    ABlast aligner(this->_maxAlignment, this->_maxPctgGap, this->_maxCtgGap);
    Contig workingCopy(ctg), origCtg(ctg);
    
    UIntType beginAlign = ctgFrame.getBegin();
    
    // alignment between pctg and ctg
    Alignment af( aligner.find_alignment((Contig)pctg, pctgPos, workingCopy, beginAlign) );
    
    beginAlign = ctg.size() - ctgFrame.getEnd();
    workingCopy = reverse_complement(origCtg);
    
    // alignment between pctg and reverse complement of ctg
    Alignment ar( aligner.find_alignment((Contig)pctg, pctgPos, workingCopy, beginAlign) );
    
    // return the best alignment found
    if(!ar.is_full())
    {
        if(!af.is_full()) return BestPctgCtgAlignment( Alignment((Contig)pctg, ctg, (size_t)pctgPos, (size_t)beginAlign), false );
        return BestPctgCtgAlignment(af,false);
    }
    else
    {
        if( !af.is_full() || ar.score() > af.score() )
        {
            ctg = workingCopy;
            return BestPctgCtgAlignment(ar,true);
        }
        
        return BestPctgCtgAlignment(af,false);
    }
}