#include <iostream>
#include <exception>
#include <boost/detail/container_fwd.hpp>

#include "pctg/PctgBuilder.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/ablast_new.hpp"
#include "alignment/banded_smith_waterman.hpp"


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


PairedContig& PctgBuilder::extendPctgWithCtgFrom(
        PairedContig& orig, 
        const Contig& ctg, 
        ContigInPctgInfo& ctgInfo, 
        const std::pair<UIntType,UIntType>& pos, 
        const std::pair<UIntType,UIntType>& gaps, 
        bool isMasterCtg) const
{
    //ctgInfo.setGaps((UIntType)gaps.first - (UIntType)gaps.second);
    
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
        //if( i < lastBases && ctg.at(ctgPos + i).base() == N ) continue;
        if( isMasterCtg || lastBases < newBases )
            orig.at(pctgPos + i) = ctg.at(ctgPos + i);
        //orig.qual(pctgPos+i) = ctg.qual(i);
    }
    
    return orig;
}


const PairedContig& PctgBuilder::extendPctgWithCtgUpto(
        PairedContig& orig, 
        const Contig& ctg, 
        const ContigInPctgInfo& ctgInfo,
        const std::pair<UIntType,UIntType>& pos,
        const UIntType pctgShift,
        bool isMasterCtg) const
{
    orig = this->shiftPctgOf(orig,pctgShift);
    
    if(isMasterCtg) orig.getMasterCtgMap()[ctgInfo.getId()] = ctgInfo;
    else orig.getSlaveCtgMap()[ctgInfo.getId()] = ctgInfo;
       
    for( UIntType i = 0; i < pos.second; i++ )
    {
        orig.at(i) = ctg.at(i);
        //orig.qual(i) = ctg.qual(i);
    }
    
    return orig;
}


PairedContig& PctgBuilder::shiftPctgOf(PairedContig& orig, const UIntType shiftSize) const
{
    orig.resize( orig.size() + shiftSize );
    
    for( IntType i = orig.size()-1; i >= shiftSize; i-- )
    {
        orig.at(i) = orig.at(i-shiftSize);
        //orig.at(shiftSize + i) = orig.at(i);
        //contigOut.qual(shiftSize + i) = orig.qual(i);
    }
    
    return shiftOf(orig,shiftSize);
}



PairedContig& PctgBuilder::addFirstContigTo(PairedContig& pctg, const IdType &ctgId) const
{
    if( pctg.size() != 0 ) throw std::invalid_argument( "Paired contig is not empty" );
    //PairedContig out(pctg.getId());
    Contig ctg = this->loadMasterContig(ctgId);
    ContigInPctgInfo ctgInfo(ctgId,ctg.size(),0);
    
    return this->extendPctgWithCtgFrom(pctg,ctg,ctgInfo,std::pair<UIntType,UIntType>(0,0),
            std::pair<UIntType,UIntType>(0,0),true);
}


PairedContig& PctgBuilder::addFirstBlockTo(PairedContig& pctg, const Block& firstBlock, const Block& lastBlock) const
{
    IdType masterCtgId = firstBlock.getMasterFrame().getContigId();
    pctg = addFirstContigTo(pctg,masterCtgId);
    
    return this->extendByBlock(pctg,firstBlock,lastBlock);
}


PairedContig& PctgBuilder::extendByBlock(PairedContig& pctg, const Block& firstBlock, const Block& lastBlock) const
{
    if( pctg.size() == 0 ) return this->addFirstBlockTo(pctg,firstBlock,lastBlock);
    
    IdType masterCtgId, slaveCtgId;
    
    masterCtgId = firstBlock.getMasterFrame().getContigId();
    slaveCtgId  = firstBlock.getSlaveFrame().getContigId();
    
    
    if( !pctg.containsMasterCtg(masterCtgId) )
    {
        if( !pctg.containsSlaveCtg(slaveCtgId) )
        {
            throw std::logic_error("Paired contig cannot be extended by this block.");
        }
        else
        {
            // merge the master contig of the block with the paired contig
            return this->mergeContig(pctg,firstBlock,lastBlock,true);
        }
    }
    else
    {
        if( !pctg.containsSlaveCtg(slaveCtgId) )
        {
            // merge the slave contig of the block with the paired contig
            return this->mergeContig(pctg,firstBlock,lastBlock,false);
        }
    }
    
    return pctg;
}


PairedContig& PctgBuilder::mergeContig(PairedContig &pctg, const Block& firstBlock, const Block& lastBlock, bool mergeMasterCtg) const
{
    IdType firstCtgId,  // identifier of the contig to be merged
           secondCtgId; // identifier of the contig already inside the paired contig
    
    if( mergeMasterCtg )
    {
        firstCtgId = firstBlock.getMasterFrame().getContigId();
        secondCtgId = firstBlock.getSlaveFrame().getContigId();
    }
    else
    {
        firstCtgId = firstBlock.getSlaveFrame().getContigId();
        secondCtgId = firstBlock.getMasterFrame().getContigId();
    }
    
    UIntType pctgPos; // position in pctg from which the alignment should be searched
    
    if( pctg.getContigInfo(secondCtgId,!mergeMasterCtg).isReversed() )
    {
        UIntType endInCtg = (mergeMasterCtg) ? 
            std::max(firstBlock.getSlaveFrame().getEnd(),lastBlock.getSlaveFrame().getEnd()) : 
            std::max(firstBlock.getMasterFrame().getEnd(),lastBlock.getMasterFrame().getEnd());
        
        pctgPos = pctg.getBasePosition(secondCtgId, endInCtg, !mergeMasterCtg);
    }
    else
    {
        UIntType beginInCtg = (mergeMasterCtg) ? 
            std::min(firstBlock.getSlaveFrame().getBegin(),lastBlock.getSlaveFrame().getBegin()) : 
            std::min(firstBlock.getMasterFrame().getBegin(),lastBlock.getMasterFrame().getBegin());
        
        pctgPos = pctg.getBasePosition(secondCtgId, beginInCtg, !mergeMasterCtg);
    }
    
    if(pctgPos >= pctg.size()) pctgPos = pctg.size()-1;
    
    // load contig that should be merged with pctg.
    Contig ctg = (mergeMasterCtg) ? this->loadMasterContig(firstCtgId) : this->loadSlaveContig(firstCtgId);    
    Frame firstFrame = (mergeMasterCtg) ? firstBlock.getMasterFrame() : firstBlock.getSlaveFrame();
    Frame lastFrame = (mergeMasterCtg) ? lastBlock.getMasterFrame() : lastBlock.getSlaveFrame();
    
    // find best alignment between ctg and pctg.
    BestPctgCtgAlignment bestAlign( this->findBestAlignment(pctg,pctg.getContigInfo(secondCtgId,!mergeMasterCtg),
                                                pctgPos, ctg, firstFrame, lastFrame) );
    
    // if the alignment is not good enough, return the paired contig
    if( bestAlign.getAlignment().homology() < MIN_HOMOLOGY ||
        bestAlign.getAlignment().length() < MIN_ALIGNMENT  ||
        bestAlign.getAlignment().score() < 0 )
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


PairedContig& PctgBuilder::mergeCtgInPos(
        PairedContig& pctg, 
        const Contig& ctg, 
        const IdType& ctgId, 
        const BestPctgCtgAlignment& bestAlign, 
        bool mergeMaster) const
{
    if( mergeMaster ) return this->mergeMasterCtgInPos(pctg,ctg,ctgId,bestAlign);
    
    return this->mergeSlaveCtgInPos(pctg,ctg,ctgId,bestAlign);
}


PairedContig& PctgBuilder::mergeMasterCtgInPos(
        PairedContig& pctg, 
        const Contig& ctg, 
        const IdType& ctgId, 
        const BestPctgCtgAlignment& bestAlign) const
{
    UIntType pctgShift = 0;
    
    //if( bestAlign.getAlignment().b_position_in_a() < 0 ) pctgShift = (-bestAlign.getAlignment().b_position_in_a());
    
    std::pair<UIntType,UIntType> pos;
    if( not first_match_pos( bestAlign.getAlignment(), pos ) ) return pctg;
    
    if( pos.second > pos.first ) pctgShift = pos.second - pos.first;
    
    ContigInPctgInfo ctgInfo( ctgId, bestAlign );
    
    if(pctgShift > 0)
    {
        ctgInfo.setPosition(0);
        pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, pos, pctgShift, true );
        pos.first += pctgShift;
    }
    
    return this->extendPctgWithCtgFrom( pctg,
                        ctg, ctgInfo, pos, std::pair<UIntType,UIntType>(0,0), true );
    
}


PairedContig& PctgBuilder::mergeSlaveCtgInPos(
        PairedContig& pctg, 
        const Contig& ctg, 
        const IdType& ctgId, 
        const BestPctgCtgAlignment& bestAlign) const
{
    UIntType pctgShift = 0; //, prevSize = pctg.size();
    
    if( bestAlign.getAlignment().b_position_in_a() < 0 )
        pctgShift = (-bestAlign.getAlignment().b_position_in_a());
    
    std::pair<UIntType,UIntType> first_pos,last_pos,gaps;
    if( not first_match_pos( bestAlign.getAlignment(), first_pos) ) return pctg;
    if( not last_match_pos( bestAlign.getAlignment(), last_pos ) ) return pctg;
    if( not gaps_before_last_match( bestAlign.getAlignment(), gaps ) ) return pctg;
    
    if( first_pos.second > first_pos.first) pctgShift = first_pos.second - first_pos.first;
    
    //pctgPos = last_pos.first + pctgShift;
    //ctgPos = last_pos.second;

    ContigInPctgInfo ctgInfo( ctgId, bestAlign );
    
    if(pctgShift > 0)
    {
        ctgInfo.setPosition(0);
        pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, first_pos, pctgShift, false );
    }
    
    last_pos.first = last_pos.first + pctgShift;
    
    //if( bestAlign.getAlignment().end_b_in_a() <= (UIntType)prevSize ) last_pos.second = ctg.size();
    
    return this->extendPctgWithCtgFrom( pctg, ctg, ctgInfo, last_pos, gaps, false );
    
}


BestPctgCtgAlignment PctgBuilder::findBestAlignment(
        const PairedContig& pctg, 
        const ContigInPctgInfo& pctgInfo, 
        const UIntType pctgPos, Contig& ctg, 
        const Frame firstFrame,
        const Frame lastFrame) const
{
    BandedSmithWaterman aligner( DEFAULT_BAND_SIZE ); //ABlast aligner; //(this->_maxAlignment, this->_maxPctgGap, this->_maxCtgGap);
    Contig workingCopy(ctg);
    
    MyAlignment af;
    UIntType beginAlignFwd = std::min( firstFrame.getBegin(), lastFrame.getBegin() );
    
    if( pctgPos < beginAlignFwd ) 
        af = aligner.find_alignment( (Contig)pctg, 0, pctg.size()-1, workingCopy, beginAlignFwd-pctgPos, workingCopy.size()-1 );
    else
        af = aligner.find_alignment( (Contig)pctg, pctgPos-beginAlignFwd, pctg.size()-1, workingCopy, 0, workingCopy.size()-1 );
    
    // alignment between pctg and ctg
    // MyAlignment af( aligner.find_alignment( (Contig)pctg, workingCopy ) );
    
    //af = aligner.find_alignment((Contig)pctg,pctgPos,pctg.size()-1,workingCopy,beginAlignFwd,workingCopy.size()-1);
    
    MyAlignment ar;
    UIntType beginAlignRev = ctg.size() - std::max( firstFrame.getEnd(), lastFrame.getEnd() );
    reverse_complement(workingCopy);
    
    if( pctgPos < beginAlignRev ) 
        ar = aligner.find_alignment( (Contig)pctg, 0, pctg.size()-1, workingCopy, beginAlignRev-pctgPos, workingCopy.size()-1 );
    else
        ar = aligner.find_alignment( (Contig)pctg, pctgPos-beginAlignRev, pctg.size()-1, workingCopy, 0, workingCopy.size()-1 );
    
    // alignment between pctg and reverse complement of ctg
    //MyAlignment ar( aligner.find_alignment((Contig)pctg, workingCopy) );
    
    //ar = aligner.find_alignment((Contig)pctg,pctgPos, pctg.size()-1,workingCopy, beginAlignRev, workingCopy.size()-1);
    
    // return the best alignment found
    if( af.homology() > ar.homology() && af.length() >= std::min(firstFrame.getLength(),lastFrame.getLength()) ) 
        return BestPctgCtgAlignment(af,false);
    
    if( ar.homology() > af.homology() && ar.length() >= std::min(firstFrame.getLength(),lastFrame.getLength()) )
    { 
        reverse_complement(ctg);
        return BestPctgCtgAlignment(ar,true); 
    }
    
    if( af.score() > ar.score() ) return BestPctgCtgAlignment( af,false );
    
    if( ar.score() > af.score() ){ reverse_complement(ctg);; return BestPctgCtgAlignment( ar,true ); }
    
    if( af.length() >= ar.length() ) return BestPctgCtgAlignment( af,false );
            
    reverse_complement(ctg);
    return BestPctgCtgAlignment( ar,true );
}