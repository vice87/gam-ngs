#include <iostream>
#include <exception>
#include <boost/detail/container_fwd.hpp>
#include <sys/stat.h>

#include "pctg/PctgBuilder.hpp"
#include "pctg/MergeInCutTailFailed.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/ablast_new.hpp"
#include "alignment/banded_smith_waterman.hpp"

//pthread_mutex_t g_badAlignMutex;
//std::ofstream g_badAlignStream;
//std::vector<uint64_t> ext_fpi;
//std::vector<uint64_t> ext_fpi_2;
//std::ofstream ext_ba_desc_stream;
//uint64_t ext_ba_count;

PctgBuilder::PctgBuilder(
        const ExtContigMemPool *masterPool,
        const ExtContigMemPool *slavePool,
        const BamTools::RefVector *masterRefVector,
        const std::vector<BamTools::RefVector> *slaveRefVector) :
                _masterPool(masterPool), _slavePool(slavePool),
                _masterRefVector(masterRefVector), _slaveRefVector(slaveRefVector),
                _maxAlignment(DEFAULT_MAX_SEARCHED_ALIGNMENT),
                _maxPctgGap(DEFAULT_MAX_GAPS),
                _maxCtgGap(DEFAULT_MAX_GAPS) {}


const Contig& PctgBuilder::loadMasterContig(const std::pair<IdType,IdType>& ctgId) const
{
    std::string refName = this->_masterRefVector->at(ctgId.second).RefName;
    return (this->_masterPool)->get( ctgId.first, refName );
}

const Contig& PctgBuilder::loadSlaveContig(const std::pair<IdType,IdType>& ctgId) const
{
    std::string refName = (this->_slaveRefVector->at(ctgId.first)).at(ctgId.second).RefName;
    return (this->_slavePool)->get( ctgId.first, refName );
}


PairedContig PctgBuilder::initByContig(const IdType& pctgId, const std::pair<IdType,IdType>& ctgId) const
{
    PairedContig pctg(pctgId);
    return this->addFirstContigTo(pctg,ctgId);
}


PairedContig& PctgBuilder::extendPctgWithCtgFrom(
        PairedContig& orig,
        const Contig& ctg,
        ContigInPctgInfo& ctgInfo,
        const std::pair<uint64_t,uint64_t>& pos,
        bool isMasterCtg) const
{
    uint64_t pctgPos = pos.first;
    uint64_t ctgPos = pos.second;

    uint64_t lastBases = (orig.size()-pctgPos);
    uint64_t newBases = (ctg.size()-ctgPos);

    if( lastBases >= newBases ) return orig;

	// resize and extend the pairedcontig
	orig.resize( pctgPos + newBases );
    for( uint64_t i = 0; i < newBases; i++ ) orig.at(pctgPos + i) = ctg.at(ctgPos + i);

    return orig;
}


PairedContig& PctgBuilder::updatePctgWithCtg(
	PairedContig& pctg,
	const Contig& ctg,
	ContigInPctgInfo& ctgInPctgInfo,
	std::pair<uint64_t,uint64_t>& start_align,
	std::pair<uint64_t,uint64_t>& end_align ) const
{
	// calcolo la differenza delle lunghezze degli intervalli dei frame coinvolti nel merging.
	uint64_t diff_a = end_align.first - start_align.first + 1;
	uint64_t diff_b = end_align.second - start_align.second + 1;
	uint64_t shift = 0;
	int64_t gap = 0;

	if( diff_a > diff_b ) // devo contrarre il pctg
	{
		shift = diff_a-diff_b;
		for( uint64_t i = end_align.first+1; i < pctg.size(); i++ ) pctg[i-shift] = pctg[i];
		pctg.resize( pctg.size()-shift );

		gap = -(int64_t(shift));
	}

	if( diff_a < diff_b ) // devo espandere il pctg
	{
		shift = diff_b-diff_a;
		pctg.resize( pctg.size()+shift );
		for( uint64_t i = pctg.size()-1; i > end_align.first+shift ; i-- ) pctg[i] = pctg[i-shift];

		gap = int64_t(shift);
	}

	// update expanded/contracted sub-sequence
	for( uint64_t i=0; i < diff_b; i++ ) pctg[start_align.first+i] = ctg[start_align.second+i];

	// update ending point of the alignment in pctg
	end_align.first = (diff_a >= diff_b) ? end_align.first-shift : end_align.first+shift;

	// update contig inside the pctg infos
	ctgInPctgInfo.addMergeGap( start_align.first, end_align.first, gap );

	return pctg;
}


PairedContig& PctgBuilder::extendPctgWithCtgUpto(
        PairedContig& orig,
        const Contig& ctg,
        const ContigInPctgInfo& ctgInfo,
        const std::pair<UIntType,UIntType>& pos,
        const UIntType pctgShift,
        bool isMasterCtg) const
{
    orig = this->shiftPctgOf(orig,pctgShift);

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



PairedContig& PctgBuilder::addFirstContigTo(PairedContig& pctg, const std::pair<IdType,IdType> &ctgId) const
{
    if( pctg.size() != 0 ) throw std::invalid_argument( "Paired contig is not empty" );
    //PairedContig out(pctg.getId());
    const Contig& ctg = this->loadMasterContig(ctgId);

	ContigInPctgInfo ctgInfo(ctgId.first, ctgId.second, ctg.size(), 0);
	pctg.getMasterCtgMap()[ctgInfo.getId()] = ctgInfo;

    return this->extendPctgWithCtgFrom(pctg,ctg,ctgInfo,std::pair<uint64_t,uint64_t>(0,0),true);
}


PairedContig& PctgBuilder::addFirstSlaveContigTo(PairedContig& pctg, const std::pair<IdType,IdType> &ctgId) const
{
	if( pctg.size() != 0 ) throw std::invalid_argument( "Paired contig is not empty" );
													   //PairedContig out(pctg.getId());
	const Contig& ctg = this->loadSlaveContig(ctgId);

	ContigInPctgInfo ctgInfo(ctgId.first, ctgId.second, ctg.size(), 0);
	pctg.getSlaveCtgMap()[ctgInfo.getId()] = ctgInfo;

	return this->extendPctgWithCtgFrom(pctg,ctg,ctgInfo,std::pair<uint64_t,uint64_t>(0,0),true);
}


PairedContig& PctgBuilder::addFirstBlockTo(PairedContig& pctg, const std::list<Block> &blocks_list, const Options &options) const
{
	const Block &firstBlock = blocks_list.front();

    std::pair<IdType,IdType> masterCtgId = std::make_pair( firstBlock.getMasterFrame().getAssemblyId(),firstBlock.getMasterFrame().getContigId() );
    pctg = addFirstContigTo(pctg,masterCtgId);

    return this->extendByBlock( pctg, blocks_list, options );
}


PairedContig& PctgBuilder::extendByBlock(PairedContig& pctg, const std::list<Block> &blocks_list, const Options &options) const
{
	if( pctg.size() == 0 ) return this->addFirstBlockTo( pctg, blocks_list, options );

    IdType masterCtgId, slaveCtgId, masterAId, slaveAId;

	const Block &block = blocks_list.front();
	const Frame &master_frame = block.getMasterFrame();
	const Frame &slave_frame = block.getSlaveFrame();

    masterAId = master_frame.getAssemblyId();
    masterCtgId = master_frame.getContigId();
    slaveAId = slave_frame.getAssemblyId();
    slaveCtgId = slave_frame.getContigId();

    if( !pctg.containsMasterCtg(masterAId,masterCtgId) )
    {
        if( !pctg.containsSlaveCtg(slaveAId,slaveCtgId) )
        {
            throw std::logic_error("Paired contig cannot be extended by this block.");
        }
        else
        {
            // merge the master contig of the block with the paired contig
            return this->mergeContig(pctg,blocks_list,true,options);
        }
    }
    else
    {
        if( !pctg.containsSlaveCtg(slaveAId,slaveCtgId) )
        {
            // merge the slave contig of the block with the paired contig
            return this->mergeContig(pctg,blocks_list,false,options);
        }
    }

    return pctg;
}


PairedContig&
PctgBuilder::mergeContig(
	PairedContig &pctg,
	const std::list<Block> &blocks_list,
	bool mergeMasterCtg,
	const Options &options) const
{
	const Block &firstBlock = blocks_list.front();
	const Block &lastBlock = blocks_list.back();

	const Frame &firstMasterFrame = firstBlock.getMasterFrame();
	const Frame &firstSlaveFrame = firstBlock.getSlaveFrame();
	const Frame &lastMasterFrame = lastBlock.getMasterFrame();
	const Frame &lastSlaveFrame = lastBlock.getSlaveFrame();

    std::pair<IdType,IdType>	firstCtgId,  // identifier of the contig to be merged
								secondCtgId; // identifier of the contig already inside the paired contig

    if( mergeMasterCtg )
    {
        firstCtgId = std::make_pair( firstMasterFrame.getAssemblyId(), firstMasterFrame.getContigId() );
        secondCtgId = std::make_pair( firstSlaveFrame.getAssemblyId(), firstSlaveFrame.getContigId() );
    }
    else
    {
        firstCtgId = std::make_pair( firstSlaveFrame.getAssemblyId(), firstSlaveFrame.getContigId() );
        secondCtgId = std::make_pair( firstMasterFrame.getAssemblyId(), firstMasterFrame.getContigId() );
    }

    int64_t startPos; 	// position in pctg from which the alignment will be searched (where the first block starts)
    int64_t endPos;		// position in pctg where the last block ends

    UIntType endInCtg = (mergeMasterCtg) ?
		std::max(firstSlaveFrame.getEnd(),lastSlaveFrame.getEnd()) :
		std::max(firstMasterFrame.getEnd(),lastMasterFrame.getEnd());

	UIntType beginInCtg = (mergeMasterCtg) ?
		std::min(firstSlaveFrame.getBegin(),lastSlaveFrame.getBegin()) :
		std::min(firstMasterFrame.getBegin(),lastMasterFrame.getBegin());

    if( pctg.getContigInfo(secondCtgId,!mergeMasterCtg).isReversed() ) // contig in pctg is reversed
    {
		startPos = pctg.getBasePosition(secondCtgId, endInCtg, !mergeMasterCtg);
		endPos = pctg.getBasePosition(secondCtgId, beginInCtg, !mergeMasterCtg);
    }
    else // contig in pctg is NOT reversed
    {
        startPos = pctg.getBasePosition(secondCtgId, beginInCtg, !mergeMasterCtg);
		endPos = pctg.getBasePosition(secondCtgId, endInCtg, !mergeMasterCtg);
    }

    //update start/end due to previous compressions/extensions during merging blocks
    int64_t gaps = 0;
    const std::list< t_merge_gap >& merge_gaps = pctg.getContigInfo(secondCtgId,!mergeMasterCtg).merge_gaps();
	for( std::list< t_merge_gap >::const_iterator mg = merge_gaps.begin(); mg != merge_gaps.end(); mg++ )
		if( mg->start < startPos ) gaps += mg->gap;

	startPos += gaps;
	endPos += gaps;

	if( startPos >= pctg.size() ) startPos = pctg.size()-1;		// for safety (it should be possible to delete this line)
	if( startPos < 0 ) startPos = 0;							// for safety (it should be possible to delete this line)
	if( endPos >= pctg.size() ) endPos = pctg.size()-1;			// for safety (it should be possible to delete this line)
	if( endPos < 0 ) endPos = 0;								// for safety (it should be possible to delete this line)

    // load contig that should be merged with pctg.
    Contig ctg = (mergeMasterCtg) ? this->loadMasterContig(firstCtgId) : this->loadSlaveContig(firstCtgId);

	//g_badAlignStream << "MERGING: " << firstMasterFrame.getContigId() << " -- " << firstSlaveFrame.getContigId() << std::endl;

    // find best alignment between ctg and pctg.
	BestPctgCtgAlignment bestAlign(
		this->findBestAlignment(pctg,pctg.getContigInfo(secondCtgId,!mergeMasterCtg),startPos,endPos,ctg,mergeMasterCtg,blocks_list)
	);

	if( bestAlign.main().homology() < MIN_HOMOLOGY )
	{
		uint64_t left_cut = pctg.getContigInfo(secondCtgId,!mergeMasterCtg).getLeftCut();
		uint64_t right_cut = pctg.getContigInfo(secondCtgId,!mergeMasterCtg).getRightCut();
		uint64_t ctg_in_pctg_size = pctg.getContigInfo(secondCtgId,!mergeMasterCtg).getSize();

		if(pctg.getContigInfo(secondCtgId,!mergeMasterCtg).isReversed()) std::swap(left_cut,right_cut);

		if( beginInCtg >= 0 && beginInCtg < left_cut && endInCtg >= 0 && endInCtg < left_cut )
			throw MergeInCutTailFailed("Merge in overwritten tail");

		if( beginInCtg >= ctg_in_pctg_size-right_cut && beginInCtg < ctg_in_pctg_size && endInCtg >= ctg_in_pctg_size-right_cut &&  endInCtg < ctg_in_pctg_size )
			throw MergeInCutTailFailed("Merge in overwritten tail");
	}

	int64_t ctg_beg = mergeMasterCtg ? std::min(firstMasterFrame.getBegin(),lastMasterFrame.getBegin()) : std::min(firstSlaveFrame.getBegin(),lastSlaveFrame.getBegin());
	int64_t ctg_end = mergeMasterCtg ? std::max(firstMasterFrame.getEnd(),lastMasterFrame.getEnd()) : std::max(firstSlaveFrame.getEnd(),lastSlaveFrame.getEnd());

	// if the alignment is not good enough, return the paired contig
    if( bestAlign.main().homology() < MIN_HOMOLOGY )
    {
		//pthread_mutex_lock(&g_badAlignMutex);
		//ext_fpi_2[ blocks_list.size() ]++;

		/*if( blocks_list.size() == 2 )
		{
			g_badAlignStream << "bad blocks alignment: ctg(0," << firstBlock.getMasterFrame().getContigId() << ") -- ctg("
			<< firstBlock.getSlaveFrame().getAssemblyId() << "," << firstBlock.getSlaveFrame().getContigId() << ")";
			g_badAlignStream << "\tstart/end in ctgInPctg = " << beginInCtg << "/" << endInCtg;
			g_badAlignStream << " start/end in pctg = " << startPos << "/" << endPos << " (size=" << pctg.size() << ")";
			g_badAlignStream << " start/end ctg = " << ctg_beg << "/" << ctg_end << " (size=" << ctg.size() << ")";
			g_badAlignStream << " gaps=" << gaps << "\n" << std::endl;

			std::stringstream ss1,ss2;
			struct stat st;

			if( !mergeMasterCtg ) ss1 << "./test3/master_" << firstBlock.getMasterFrame().getContigId() << "_" << ext_ba_count;
			else ss1 << "./test3/slave_" << firstBlock.getSlaveFrame().getContigId() << "_" << ext_ba_count;

			if( mergeMasterCtg ) ss2 << "./test3/master_" << firstBlock.getMasterFrame().getContigId() << "_" << ext_ba_count;
			else ss2 << "./test3/slave_" << firstBlock.getSlaveFrame().getContigId() << "_" << ext_ba_count;

			if( stat(ss1.str().c_str(),&st) == 0 || stat(ss2.str().c_str(),&st) == 0 )
			{
				ext_ba_count++;

				ss1.str( std::string("") );
				ss1.clear();
				ss2.str( std::string("") );
				ss2.clear();

				if( !mergeMasterCtg ) ss1 << "./test3/master_" << firstBlock.getMasterFrame().getContigId() << "_" << ext_ba_count;
				else ss1 << "./test3/slave_" << firstBlock.getSlaveFrame().getContigId() << "_" << ext_ba_count;

				if( mergeMasterCtg ) ss2 << "./test3/master_" << firstBlock.getMasterFrame().getContigId() << "_" << ext_ba_count;
				else ss2 << "./test3/slave_" << firstBlock.getSlaveFrame().getContigId() << "_" << ext_ba_count;
			}

			std::ofstream ba_out( ss1.str().c_str() );

			if( !mergeMasterCtg ) ba_out << ">master_" << firstBlock.getMasterFrame().getContigId() << std::endl;
			else ba_out << ">slave_" << firstBlock.getSlaveFrame().getContigId() << std::endl;

			for( int64_t i = (startPos>=200 ? startPos-200 : 0); i <= endPos + 200; i++ )
			{
				if( i < pctg.size() ) ba_out << char(pctg[i]);
			}
			ba_out << std::endl;
			ba_out.close();

			ba_out.open( ss2.str().c_str() );

			if( mergeMasterCtg ) ba_out << ">master_" << firstBlock.getMasterFrame().getContigId() << std::endl;
			else ba_out << ">slave_" << firstBlock.getSlaveFrame().getContigId() << std::endl;

			for( int64_t i = (ctg_beg >= 200 ? ctg_beg-200 : 0); i <= ctg_end+ 200; i++ )
			{
				if( i < ctg.size() ) ba_out << char(ctg[i]);
			}
			ba_out << std::endl;
			ba_out.close();

			ext_ba_desc_stream << "master_" << firstBlock.getMasterFrame().getContigId() << "_" << ext_ba_count << std::endl;
			ext_ba_desc_stream << "slave_" << firstBlock.getSlaveFrame().getContigId() << "_" << ext_ba_count << std::endl;
		}*/

		//pthread_mutex_unlock(&g_badAlignMutex);
        return pctg;
    }
    else if( bestAlign.left().homology() < MIN_HOMOLOGY_2 || bestAlign.right().homology() < MIN_HOMOLOGY_2 )
	{
		//pthread_mutex_lock(&g_badAlignMutex);
		//ext_fpi_2[ blocks_list.size() ]++;

		//g_badAlignStream << "good blocks alignment:\t(0," << firstBlock.getMasterFrame().getContigId() << ")\t("
		//	<< firstBlock.getSlaveFrame().getAssemblyId() << "," << firstBlock.getSlaveFrame().getContigId() << ")"
		//	<< "\tidentity=" << bestAlign.main().homology() << std::endl;

		//g_badAlignStream << "bad tail alignment:";

		//if( bestAlign.left().homology() < MIN_HOMOLOGY_2 ) g_badAlignStream << "\tLEFT(" << bestAlign.left().homology() << ")";
		//if( bestAlign.right().homology() < MIN_HOMOLOGY_2 ) g_badAlignStream << "\tRIGHT(" << bestAlign.right().homology() << ")";

		//g_badAlignStream << "\n" << std::endl;

		//pthread_mutex_unlock(&g_badAlignMutex);
		return pctg;
	}

	//pthread_mutex_lock(&g_badAlignMutex);
	//ext_fpi[ blocks_list.size() ]++;
	//pthread_mutex_unlock(&g_badAlignMutex);

    // avoid to overlap primary assembly contigs
    /*if( mergeMasterCtg )
    {
        if( sameAssemblyCtgsOverlapedBy(pctg,ctg,bestAlign.main().b_position_in_a(), mergeMasterCtg) )
        {
            //std::cout << "No overlap for primary assembly constraint disattended." << std::endl << std::flush;
            return pctg;
            throw ConstraintsDisattended("No overlap for primary assembly constraint disattended.");
        }
    }*/

    // merge ctg inside pctg, using alignment informations.
    return this->mergeCtgInPos(pctg,ctg,secondCtgId,firstCtgId,bestAlign,mergeMasterCtg);
}


PairedContig& PctgBuilder::mergeCtgInPos(
        PairedContig& pctg,
        const Contig& ctg,
		const std::pair<IdType,IdType>& ctgInPctgId,
        const std::pair<IdType,IdType>& ctgId,
        const BestPctgCtgAlignment& alignment,
        bool mergeMaster) const
{
	// following: previous code - DON'T DELETE
	/*if( mergeMaster ) return this->mergeMasterCtgInPos(pctg,ctg,ctgId,bestAlign);
	return this->mergeSlaveCtgInPos(pctg,ctg,ctgId,bestAlign);*/

	// new way of merging
	uint64_t pctg_shift = 0;
	std::pair<uint64_t,uint64_t> start_align, end_align;
	ContigInPctgInfo ctgInfo( ctgId.first, ctgId.second, alignment );

	ContigInPctgInfo& ctgInPctgInfo = mergeMaster ? pctg.getSlaveCtgMap()[ctgInPctgId] : pctg.getMasterCtgMap()[ctgInPctgId];

	first_match_pos( alignment.main(), start_align );
	last_match_pos( alignment.main(), end_align );

	// if possible extend the pctg from the left
	if( start_align.second > start_align.first )
	{
		pctg_shift = start_align.second - start_align.first;
		if(start_align.first > 0) ctgInPctgInfo.setLeftCut(start_align.first);
	}
	else
	{
		if(start_align.second > 0) ctgInfo.setLeftCut(start_align.second);
	}

	if( pctg_shift > 0 )
	{
		ctgInfo.setPosition(0);
		pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, start_align, pctg_shift, mergeMaster );
		start_align.first += pctg_shift;
		end_align.first += pctg_shift;
	}

	// if ctg contains less Ns, update the central part of pctg
	uint64_t pctg_n = 0, ctg_n = 0;
	for( uint64_t i = start_align.first; i <= end_align.first; i++ ) if( char(pctg[i]) == 'N' ) pctg_n++;
	for( uint64_t i = start_align.second; i <= end_align.second; i++ ) if( char(ctg[i]) == 'N' ) ctg_n++;

	if( ctg_n <= pctg_n )
	{
		pctg = this->updatePctgWithCtg( pctg, ctg, ctgInPctgInfo, start_align, end_align );

	}
	else // update ctgInfo due to compression/extension
	{
		uint64_t diff_a = end_align.first - start_align.first + 1;
		uint64_t diff_b = end_align.second - start_align.second + 1;
		int64_t gap = diff_a - diff_b;

		ctgInfo.addMergeGap( start_align.second, end_align.second, gap );
	}

	if( pctg.size()-end_align.first >= ctg.size()-end_align.second ) // ctg right tail cut
	{
		ctgInPctgInfo.setRightCut(ctg.size()-end_align.second-1);
	}
	else // ctg in pctg tail cut
	{
		ctgInPctgInfo.setRightCut(pctg.size()-end_align.first-1);
	}

	pctg = this->extendPctgWithCtgFrom( pctg, ctg, ctgInfo, end_align, mergeMaster);

	// update contig info in pctg
	if(mergeMaster) pctg.getMasterCtgMap()[ctgInfo.getId()] = ctgInfo;
		else pctg.getSlaveCtgMap()[ctgInfo.getId()] = ctgInfo;

	return pctg;
}


/*PairedContig& PctgBuilder::mergeMasterCtgInPos(
        PairedContig& pctg,
        const Contig& ctg,
        const std::pair<IdType,IdType>& ctgId,
        const BestPctgCtgAlignment& bestAlign) const
{
    uint64_t pctgShift = 0;

    //if( bestAlign.getAlignment().b_position_in_a() < 0 ) pctgShift = (-bestAlign.getAlignment().b_position_in_a());

    std::pair<UIntType,UIntType> pos;
    if( not first_match_pos( bestAlign.main(), pos ) ) return pctg;

    if( pos.second > pos.first ) pctgShift = pos.second - pos.first;

    ContigInPctgInfo ctgInfo( ctgId.first, ctgId.second, bestAlign );

    if(pctgShift > 0)
    {
        ctgInfo.setPosition(0);
        pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, pos, pctgShift, true );
        pos.first += pctgShift;
    }

    return this->extendPctgWithCtgFrom( pctg, ctg, ctgInfo, pos, std::pair<UIntType,UIntType>(0,0), true );

}


PairedContig& PctgBuilder::mergeSlaveCtgInPos(
        PairedContig& pctg,
        const Contig& ctg,
        const std::pair<IdType,IdType>& ctgId,
        const BestPctgCtgAlignment& bestAlign) const
{
    UIntType pctgShift = 0; //, prevSize = pctg.size();

    if( bestAlign.main().b_position_in_a() < 0 )
        pctgShift = (-bestAlign.main().b_position_in_a());

    std::pair<UIntType,UIntType> first_pos,last_pos,gaps;
    if( not first_match_pos( bestAlign.main(), first_pos) ) return pctg;
    if( not last_match_pos( bestAlign.main(), last_pos ) ) return pctg;
    if( not gaps_before_last_match( bestAlign.main(), gaps ) ) return pctg;

    if( first_pos.second > first_pos.first) pctgShift = first_pos.second - first_pos.first;

    //pctgPos = last_pos.first + pctgShift;
    //ctgPos = last_pos.second;

    ContigInPctgInfo ctgInfo( ctgId.first, ctgId.second, bestAlign );

    if(pctgShift > 0)
    {
        ctgInfo.setPosition(0);
        pctg = this->extendPctgWithCtgUpto( pctg, ctg, ctgInfo, first_pos, pctgShift, false );
    }

    last_pos.first = last_pos.first + pctgShift;

    //if( bestAlign.getAlignment().end_b_in_a() <= (UIntType)prevSize ) last_pos.second = ctg.size();

    return this->extendPctgWithCtgFrom( pctg, ctg, ctgInfo, last_pos, gaps, false );

}*/


BestPctgCtgAlignment PctgBuilder::findBestAlignment(
        const PairedContig &pctg,
		const ContigInPctgInfo& pctgInfo,
        const uint64_t startPos,
		const uint64_t endPos,
		Contig &ctg,
		bool isMasterCtg,
		const std::list<Block> &blocks_list) const
{
	bool reversed = pctgInfo.isReversed(); // whether the contig inside the pctg is reversed or not.
	uint64_t con_evid = 0, dis_evid = 0;
	uint64_t mf_len = 0, sf_len = 0; // sum of master/slave frames lengths

	// compute the probability of ctg to be reverse complemented respect to the contig inside the pctg
	for( std::list<Block>::const_iterator b = blocks_list.begin(); b != blocks_list.end(); b++ )
	{
		const Frame& mf = b->getMasterFrame();
		const Frame& sf = b->getSlaveFrame();

		mf_len += mf.getLength();
		sf_len += sf.getLength();

		if( mf.getStrand() != sf.getStrand() ) dis_evid += b->getReadsNumber(); // discordant frames
			else con_evid += b->getReadsNumber(); // concordant frames
	}

	uint64_t min_align_len = std::min(mf_len,sf_len);
	double con_prob = double(con_evid) / double(con_evid+dis_evid);

	//g_badAlignStream << "evidence = " << con_prob << " reversed = " << (reversed ? "true" : "false")  << std::endl;

	// first & last blocks references
	const Block &first_block = blocks_list.front();
	const Block &last_block = blocks_list.back();

	// first & last frames (of ctg) references
	const Frame &ctgFirstFrame = (isMasterCtg) ? first_block.getMasterFrame() : first_block.getSlaveFrame();
	const Frame &ctgLastFrame = (isMasterCtg) ? last_block.getMasterFrame() : last_block.getSlaveFrame();

	// blocks interval in the contig to be merged
	uint64_t startInCtg = std::min( ctgFirstFrame.getBegin(), ctgLastFrame.getBegin() );
	uint64_t endInCtg = std::max( ctgFirstFrame.getEnd(), ctgLastFrame.getEnd() );
	uint64_t tempInCtg;

	// compute tails / min alignment length thresholds
	const Frame& mf = first_block.getMasterFrame();
	const Frame& sf = first_block.getSlaveFrame();
	int32_t m_len = 0.3 * (_masterRefVector->at(mf.getContigId())).RefLength;
	int32_t s_len = 0.3 * (_slaveRefVector->at(sf.getAssemblyId())).at(sf.getContigId()).RefLength;
	int32_t threshold = std::min( 500, std::min( m_len, s_len ) );

	uint64_t pctg_int = endPos >= startPos ? endPos-startPos : 0;
	uint64_t ctg_int = endInCtg >= startInCtg ? endInCtg-startInCtg : 0;
	uint64_t pctg_ctg_int_diff = pctg_int >= ctg_int ? pctg_int-ctg_int : ctg_int-pctg_int;
	pctg_ctg_int_diff += 0.2 * pctg_ctg_int_diff;

	uint64_t bsw_band = (pctg_ctg_int_diff > 1000) ? DEFAULT_BAND_SIZE : std::max( uint64_t(DEFAULT_BAND_SIZE), pctg_ctg_int_diff );

	BandedSmithWaterman aligner( bsw_band ); //ABlast aligner; //(this->_maxAlignment, this->_maxPctgGap, this->_maxCtgGap);
	MyAlignment align, bad_align(0), good_align(100);

	bool good_align_found = false;
	bool is_ctg_rev = false;

	// contigs more likely have the same orientation
	if( (con_prob >= 0.5 && !reversed) || (con_prob <= 0.5 && reversed) )
	{
		align = aligner.find_alignment( (Contig)pctg, startPos, pctg.size()-1, ctg, startInCtg, endInCtg );

		if( this->is_good( align, threshold ) )
		{
			good_align_found = true;
			is_ctg_rev = false;
		}
		else
		{
			/*pthread_mutex_lock(&g_badAlignMutex);
			g_badAlignStream << mf.getContigId() << "/" << sf.getContigId() << " hom=" << align.homology() << " len=" << align.length() << "/" << threshold;
			g_badAlignStream << " pctg_pos=" << startPos << " ctg_pos=" << startInCtg << "/" << endInCtg << std::endl;
			pthread_mutex_unlock(&g_badAlignMutex);*/

			// else, try reversing the contig
			reverse_complement(ctg);

			// update start and end positions of the blocks
			tempInCtg = startInCtg;
			startInCtg = ctg.size() - endInCtg - 1;
			endInCtg = ctg.size() - tempInCtg - 1;

			align = aligner.find_alignment( (Contig)pctg, startPos, pctg.size()-1, ctg, startInCtg, endInCtg );

			// if align is good return it
			if( this->is_good( align, threshold ) ){ good_align_found = true; is_ctg_rev = true; }
			/*else
			{
				pthread_mutex_lock(&g_badAlignMutex);
				g_badAlignStream << mf.getContigId() << "/" << sf.getContigId() << " hom=" << align.homology() << " len=" << align.length() << "/" << threshold;
				g_badAlignStream << " pctg_pos=" << startPos << " ctg_pos=" << startInCtg << "/" << endInCtg << std::endl;
				pthread_mutex_unlock(&g_badAlignMutex);
			}*/
		}
	}

	// contigs more likely have opposite orientations
	if( (con_prob > 0.5 && reversed) || (con_prob < 0.5 && !reversed))
	{
		reverse_complement(ctg);

		// update start and end positions of the blocks
		tempInCtg = startInCtg;
		startInCtg = ctg.size() - endInCtg - 1;
		endInCtg = ctg.size() - tempInCtg - 1;

		align = aligner.find_alignment( (Contig)pctg, startPos, pctg.size()-1, ctg, startInCtg, endInCtg );

		if( this->is_good( align, threshold ) )
		{
			good_align_found = true;
			is_ctg_rev = true;
		}
		else
		{
			/*pthread_mutex_lock(&g_badAlignMutex);
			g_badAlignStream << mf.getContigId() << "/" << sf.getContigId() << " hom=" << align.homology() << " len=" << align.length() << "/" << threshold;
			g_badAlignStream << " pctg_pos=" << startPos << " ctg_pos=" << startInCtg << "/" << endInCtg << std::endl;
			pthread_mutex_unlock(&g_badAlignMutex);*/

			// restore original contig
			reverse_complement(ctg);

			// update start and end positions of the blocks
			tempInCtg = startInCtg;
			startInCtg = ctg.size() - endInCtg - 1;
			endInCtg = ctg.size() - tempInCtg - 1;

			align = aligner.find_alignment( (Contig)pctg, startPos, pctg.size()-1, ctg, startInCtg, endInCtg );

			if( this->is_good( align, threshold ) ){ good_align_found = true; is_ctg_rev = false; }
			/*else
			{
				pthread_mutex_lock(&g_badAlignMutex);
				g_badAlignStream << mf.getContigId() << "/" << sf.getContigId() << " hom=" << align.homology() << " len=" << align.length() << "/" << threshold;
				g_badAlignStream << " pctg_pos=" << startPos << " ctg_pos=" << startInCtg << "/" << endInCtg << std::endl;
				pthread_mutex_unlock(&g_badAlignMutex);
			}*/
		}
	}

	// if the alignments computed were all bad, return a bad alignment to interrupt the merging
	if( !good_align_found ) return BestPctgCtgAlignment(bad_align,is_ctg_rev);

	// retrieve the starting/ending points of the alignment
	std::pair<uint64_t,uint64_t> start_align, end_align;
	first_match_pos( align, start_align );
	last_match_pos( align, end_align );

	// compute pctg tails (i1,i2) and ctg tails (j1,j2)
	uint64_t i1 = start_align.first;
	uint64_t i2 = pctg.size() - end_align.first - 1;
	uint64_t j1 = start_align.second;
	uint64_t j2 = ctg.size() - end_align.second - 1;

	if( std::min(i1,j1) < threshold && std::min(i2,j2) < threshold ) return BestPctgCtgAlignment(align,is_ctg_rev);

	MyAlignment left(100),right(100);
	bool left_rev, right_rev;

	if( std::min(i1,j1) >= threshold ) // left tail alignment
	{
		if( i1 < j1 ) // pctg left tail < ctg left tail
		{
			left = aligner.find_alignment( ctg, start_align.second-start_align.first, start_align.second-1, (Contig)pctg, 0, start_align.first-1, false, true );
			left_rev = true;
		}
		else	// pctg left tail >= ctg left tail
		{
			left = aligner.find_alignment( (Contig)pctg, start_align.first-start_align.second, start_align.first-1, ctg, 0, start_align.second-1, false, true );
			left_rev = false;
		}
	}

	if( std::min(i2,j2) >= threshold ) // right tail alignment
	{
		if( i2 < j2 ) // pctg right tail < ctg right tail
		{
			Contig right_tail = chop_begin( ctg, end_align.second+1 );
			right = aligner.find_alignment( right_tail, 0, right_tail.size()-1, (Contig)pctg, end_align.first+1, pctg.size()-1, true, false );
			right_rev = true;
		}
		else	// pctg right tail >= ctg right tail
		{
			Contig right_tail = chop_begin( (Contig)pctg, end_align.first+1 );
			right = aligner.find_alignment( right_tail, 0, right_tail.size()-1, ctg, end_align.second+1, ctg.size()-1, true, false );
			right_rev = false;
		}
	}

	return BestPctgCtgAlignment(align,is_ctg_rev,left,right,left_rev,right_rev);
}


bool PctgBuilder::is_good( const MyAlignment &align, uint64_t min_align_len ) const
{
	return (align.homology() >= MIN_HOMOLOGY && align.length() >= min_align_len); // && align.score() > 0
}

