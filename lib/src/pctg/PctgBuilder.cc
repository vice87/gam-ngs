#include <iostream>
#include <exception>
#include <boost/detail/container_fwd.hpp>
#include <sys/stat.h>

#include "pctg/PctgBuilder.hpp"
#include "pctg/MergeInCutTailFailed.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/ablast.hpp"
#include "alignment/banded_smith_waterman.hpp"

std::ofstream g_badAlignStream;
//std::ofstream ext_ba_desc_stream;
//uint64_t ext_ba_count;

extern MultiBamReader masterBam;
extern MultiBamReader masterMpBam;
extern MultiBamReader slaveBam;
extern MultiBamReader slaveMpBam;

pthread_mutex_t g_badAlignMutex;
std::vector<uint64_t> ext_fpi;
std::vector<uint64_t> ext_fpi_2;

PctgBuilder::PctgBuilder(
	ThreadedBuildPctg *tbp,
	const RefSequence *masterRef,
	const RefSequence *slaveRef) :
		_tbp(tbp),
		_masterRef(masterRef),
		_slaveRef(slaveRef),
		_maxAlignment(DEFAULT_MAX_SEARCHED_ALIGNMENT),
		_maxPctgGap(DEFAULT_MAX_GAPS),
		_maxCtgGap(DEFAULT_MAX_GAPS)
{}


const Contig& PctgBuilder::loadMasterContig(const int32_t ctgId) const
{
    return *(this->_masterRef->at(ctgId)).Sequence;

	//std::string refName = this->_masterRefVector->at(ctgId.second).RefName;
    //return (this->_masterPool)->get( ctgId.first, refName );
}

const Contig& PctgBuilder::loadSlaveContig(const int32_t ctgId) const
{
    return *(this->_slaveRef->at(ctgId)).Sequence;

	//std::string refName = (this->_slaveRefVector->at(ctgId.first)).at(ctgId.second).RefName;
    //return (this->_slavePool)->get( ctgId.first, refName );
}


PairedContig PctgBuilder::initByContig(const IdType& pctgId, const int32_t ctgId) const
{
    PairedContig pctg(pctgId);
    return this->addFirstContigTo(pctg,ctgId);
}


void PctgBuilder::appendMasterToPctg( PairedContig &pctg, MergeStruct &ctgInfo, int64_t start, int64_t end )
{
	int32_t id = ctgInfo.id;
	Contig *ctg = ctgInfo.ctg;
	bool reversed = ctgInfo.reversed;

	if( end < start || start < 0 || end >= ctg->size() ) return;

	pctg.addMasterCtgId( ctgInfo.id );

	int64_t idx = pctg.size();
	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg->at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_back( CtgInPctgInfo( id, start, end, reversed, true ) );
}


void PctgBuilder::appendSlaveToPctg( PairedContig &pctg, MergeStruct &ctgInfo, int64_t start, int64_t end )
{
	int32_t id = ctgInfo.id;
	Contig *ctg = ctgInfo.ctg;
	bool reversed = ctgInfo.reversed;

	if( end < start || start < 0 || end >= ctg->size() ) return;

	pctg.addSlaveCtgId( ctgInfo.id );

	int64_t idx = pctg.size();
	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg->at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_back( CtgInPctgInfo( id, start, end, reversed, false ) );
}


void PctgBuilder::prependMasterToPctg( PairedContig &pctg, MergeStruct &ctgInfo, int64_t start, int64_t end )
{
	int32_t id = ctgInfo.id;
	Contig *ctg = ctgInfo.ctg;
	bool reversed = ctgInfo.reversed;

	if( end < start || start < 0 || end >= ctg->size() ) return;

	pctg.addMasterCtgId( ctgInfo.id );

	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = pctg.size()-1; i >= newBases; i-- ) pctg.at(i) = pctg.at(i-newBases);

	int64_t idx = 0;
	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg->at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_front( CtgInPctgInfo( id, start, end, reversed, true ) );
}

void PctgBuilder::prependSlaveToPctg( PairedContig &pctg, MergeStruct &ctgInfo, int64_t start, int64_t end )
{
	int32_t id = ctgInfo.id;
	Contig *ctg = ctgInfo.ctg;
	bool reversed = ctgInfo.reversed;

	if( end < start || start < 0 || end >= ctg->size() ) return;

	pctg.addSlaveCtgId( ctgInfo.id );

	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = pctg.size()-1; i >= newBases; i-- ) pctg.at(i) = pctg.at(i-newBases);

	int64_t idx = 0;
	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg->at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_front( CtgInPctgInfo( id, start, end, reversed, false ) );
}


void PctgBuilder::appendBlocksRegionToPctg( PairedContig &pctg,
											 MergeStruct &masterInfo, int64_t masterStart, int64_t masterEnd,
											 MergeStruct &slaveInfo, int64_t slaveStart, int64_t slaveEnd )
{
	pctg.addMasterCtgId( masterInfo.id );
	pctg.addSlaveCtgId( slaveInfo.id );

	int64_t master_int = (masterEnd >= masterStart) ? masterEnd - masterStart + 1 : 0;
	int64_t slave_int = (slaveEnd >= slaveStart) ? slaveEnd - slaveStart + 1 : 0;

	int64_t large_int = (master_int >= slave_int) ? master_int : slave_int;
	int64_t small_int = (master_int >= slave_int) ? slave_int : master_int;

	// if regions are similar in length append master's blocks region and return
	if( small_int >= 0.97*large_int ) return appendMasterToPctg( pctg, masterInfo, masterStart, masterEnd );

	std::vector<double> masterScore, slaveScore;

	masterScore = computeZScore( masterBam, masterInfo.id, masterStart, masterEnd );
	slaveScore = computeZScore( slaveBam, slaveInfo.id, slaveStart, slaveEnd );

	size_t masterEvid = 0;
	size_t slaveEvid = 0;

	for( size_t i=0; i < masterScore.size(); i++ )
	{
		double m_score = masterScore[i]; if( m_score < 0 ) m_score = -m_score;
		double s_score =  slaveScore[i]; if( s_score < 0 ) s_score = -s_score;

		if( s_score < m_score && s_score != 0 ) slaveEvid++; else if( s_score < m_score ) masterEvid++;
		if( m_score < s_score && m_score != 0 ) masterEvid++; else if( m_score < s_score ) slaveEvid++;
	}

	if( masterEvid >= slaveEvid ) return appendMasterToPctg( pctg, masterInfo, masterStart, masterEnd );

	return appendSlaveToPctg( pctg, slaveInfo, slaveStart, slaveEnd );
}


void PctgBuilder::prependBlocksRegionToPctg( PairedContig &pctg,
											 MergeStruct &masterInfo, int64_t masterStart, int64_t masterEnd,
											 MergeStruct &slaveInfo, int64_t slaveStart, int64_t slaveEnd )
{
	pctg.addMasterCtgId( masterInfo.id );
	pctg.addSlaveCtgId( slaveInfo.id );

	int64_t master_int = (masterEnd >= masterStart) ? masterEnd - masterStart + 1 : 0;
	int64_t slave_int = (slaveEnd >= slaveStart) ? slaveEnd - slaveStart + 1 : 0;

	int64_t large_int = (master_int >= slave_int) ? master_int : slave_int;
	int64_t small_int = (master_int >= slave_int) ? slave_int : master_int;

	// if regions are similar in length append master's blocks region and return
	if( small_int >= 0.97*large_int ) return prependMasterToPctg( pctg, masterInfo, masterStart, masterEnd );

	std::vector<double> masterScore, slaveScore;

	masterScore = computeZScore( masterBam, masterInfo.id, masterStart, masterEnd );
	slaveScore = computeZScore( slaveBam, slaveInfo.id, slaveStart, slaveEnd );

	size_t masterEvid = 0;
	size_t slaveEvid = 0;

	for( size_t i=0; i < masterScore.size(); i++ )
	{
		double m_score = masterScore[i]; if( m_score < 0 ) m_score = -m_score;
		double s_score =  slaveScore[i]; if( s_score < 0 ) s_score = -s_score;

		//if( m_score <= s_score ) masterEvid++; else slaveEvid++;
		if( s_score < m_score && s_score != 0 ) slaveEvid++; else if( s_score < m_score ) masterEvid++;
		if( m_score < s_score && m_score != 0 ) masterEvid++; else if( m_score < s_score ) slaveEvid++;
	}

	if( masterEvid >= slaveEvid ) return prependMasterToPctg( pctg, masterInfo, masterStart, masterEnd );

	return prependSlaveToPctg( pctg, slaveInfo, slaveStart, slaveEnd );
}



void PctgBuilder::reverseCoordinates( std::pair<int64_t,int64_t> &coords, size_t ctg_size )
{
	int64_t tempPos = coords.first;
	coords.first = ctg_size - coords.second - 1;
	coords.second = ctg_size - tempPos - 1;
}


void PctgBuilder::mergeContigs( std::list<PairedContig> &pctgList, MergeDescriptorLists &mergeLists )
{
	for( MergeDescriptorLists::iterator it = mergeLists.begin(); it != mergeLists.end(); ++it )
	{
		if( it->size() == 0 ) continue;

		std::list< MergeDescriptor >::iterator last = --(it->end());

		if( last->mergeType == FORK_MERGE && last->forkType == UNKNOWN )
		{
			int32_t mId = last->masterId;
			std::list< MergeDescriptor >::iterator md = last;
			std::list< MergeDescriptor >::iterator del = last;

			while( md != it->begin() )
			{
				--md;
				if( del->masterId == md->masterId ) break;
				del = md;
			}

			it->erase(del,it->end());
		}

		if( it->size() == 0 ) continue;

		this->mergeContigs( pctgList, *it );
	}
}


void PctgBuilder::mergeContigs( std::list<PairedContig> &pctgList, std::list<MergeDescriptor> &mergeList )
{
	if( mergeList.size() == 0 ) return;

	PairedContig pctg;

	bool fwdMerge = true;
	bool prevFwdMerge = true;
	bool curAlignFailed = false;
	bool prevAlignFailed = false;
	MergeStruct prevMaster, prevSlave, curMaster, curSlave;

	prevMaster.id = -1;
	prevSlave.id  = -1;
	//curMaster.id = -1;
	//curSlave.id = -1;

	std::list<MergeDescriptor>::iterator it, cur, next;

	it = mergeList.begin();
	while( it != mergeList.end() )
	{
		cur = it;
		next = ++it;

		if( cur == mergeList.begin() ) // first merge
		{
			curMaster.id = cur->masterId;
			curMaster.ctg = new Contig( this->loadMasterContig( cur->masterId ) );
			curMaster.reversed = false;
			curMaster.pos = 0;

			curSlave.id = cur->slaveId;
			curSlave.ctg = new Contig( this->loadSlaveContig( cur->slaveId ) );
			curSlave.reversed = (cur->align)->isCtgReversed();
			curSlave.pos = 0;

			size_t masterSize = (curMaster.ctg)->size();
			size_t slaveSize = (curSlave.ctg)->size();

			if( curSlave.reversed ) // if slave reversed, recompute coordinates
			{
				reverse_complement( *(curSlave.ctg) );

				this->reverseCoordinates( cur->slaveMainAlign, slaveSize );
				this->reverseCoordinates( cur->slavePos, slaveSize );
				this->reverseCoordinates( cur->slaveBlocks, slaveSize );
			}

			curMaster.limits = cur->masterPos;
			curSlave.limits = cur->slavePos;

			if( next != mergeList.end() ) // not last merge
			{
				std::pair<int64_t,int64_t> nextSlave = next->slaveBlocks;

				if( curSlave.reversed ) this->reverseCoordinates( nextSlave, slaveSize );

				if( cur->masterId == next->masterId )
				{
					fwdMerge = (cur->masterBlocks).first <= (next->masterBlocks).first;
					prevFwdMerge = fwdMerge;
				}
				else // cur->slaveId == next->slaveId
				{
					if( !curSlave.reversed ) fwdMerge = (cur->slaveBlocks).first <= nextSlave.first;
					if( curSlave.reversed ) fwdMerge = (cur->slaveBlocks).second <= nextSlave.second;
					prevFwdMerge = fwdMerge;
				}
			}

			curAlignFailed = this->isAlignmentFailed( *(cur->align), cur->masterPos, cur->slavePos, curMaster.reversed );

			if( !curAlignFailed ) // good alignment
			{
				// estendi il pctg in base alla coda più lunga fino alla regione dei blocchi (compresa)

				if( fwdMerge ) // merge a partire "dall'inizio" del pctg
				{
					// la coda più lunga è quella del master
					if( cur->masterMainAlign.first - cur->masterPos.first >= cur->slaveMainAlign.first - cur->slavePos.first )
					{
						this->appendMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterMainAlign.first - 1 );
						pctg.addMasterCtgId( curMaster.id );
					}
					else // la coda più lunga è dello slave
					{
						this->appendSlaveToPctg( pctg, curSlave, cur->slavePos.first, cur->slaveMainAlign.first - 1 );
						pctg.addSlaveCtgId( curSlave.id );
					}

					this->appendBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
													curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

					curMaster.pos = cur->masterMainAlign.second + 1;
					curSlave.pos = cur->slaveMainAlign.second + 1;
				}
				else // merge a partire dalla "fine" del pctg
				{
					if( cur->masterPos.second - cur->masterMainAlign.second >= cur->slavePos.second - cur->slaveMainAlign.second )
					{
						this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.first+1, cur->masterPos.second );
						pctg.addMasterCtgId( curMaster.id );
					}
					else
					{
						this->prependSlaveToPctg( pctg, curSlave, cur->slaveMainAlign.second+1, cur->slavePos.second );
						pctg.addSlaveCtgId( curSlave.id );
					}

					this->prependBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
													curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

					curMaster.pos = cur->masterMainAlign.first - 1;
					curSlave.pos = cur->slaveMainAlign.first - 1;
				}

				//prevAlignFailed = false;
			}
			else // bad alignment
			{
				if( fwdMerge )
				{
					this->appendMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterPos.second );
					pctg.addMasterCtgId( curMaster.id );
					curMaster.pos = cur->masterPos.first; //cur->masterPos.second + 1;
				}
				else
				{
					this->prependMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterPos.second );
					pctg.addMasterCtgId( curMaster.id );
					curMaster.pos = cur->masterPos.second; //cur->masterPos.first - 1;
				}

				//prevAlignFailed = true;
			}
		}
		else // not first merge
		{
			curMaster.id = cur->masterId;
			curSlave.id = cur->slaveId;

			if( prevMaster.id == curMaster.id )
			{
				curMaster.ctg = prevMaster.ctg;
				curSlave.ctg = new Contig( this->loadSlaveContig( cur->slaveId ) );
				curSlave.reversed = (prevMaster.reversed && !(cur->align)->isCtgReversed() || !prevMaster.reversed && (cur->align)->isCtgReversed() );
				curSlave.pos = 0;

				if( curSlave.reversed ) reverse_complement( *(curSlave.ctg) );
			}
			else // prevSlave.id == curSlave.id
			{
				curSlave = prevSlave;
				curMaster.ctg = new Contig( this->loadMasterContig( cur->masterId ) );
				curMaster.reversed = (prevSlave.reversed && !(cur->align)->isCtgReversed() || !prevSlave.reversed && (cur->align)->isCtgReversed() );
				curMaster.pos = 0;

				if( curMaster.reversed ) reverse_complement( *(curMaster.ctg) );
			}

			size_t masterSize = (curMaster.ctg)->size();
			size_t slaveSize = (curSlave.ctg)->size();

			if( curMaster.reversed ) // if master reversed, recompute coordinates
			{
				this->reverseCoordinates( cur->masterMainAlign, masterSize );
				this->reverseCoordinates( cur->masterPos, masterSize );
				this->reverseCoordinates( cur->masterBlocks, masterSize );
			}

			if( curSlave.reversed ) // if slave reversed, recompute coordinates
			{
				this->reverseCoordinates( cur->slaveMainAlign, slaveSize );
				this->reverseCoordinates( cur->slavePos, slaveSize );
				this->reverseCoordinates( cur->slaveBlocks, slaveSize );
			}

			curAlignFailed = this->isAlignmentFailed( *(cur->align), cur->masterPos, cur->slavePos, curMaster.reversed );

			curMaster.limits = cur->masterPos;
			curSlave.limits = cur->slavePos;

			if( next != mergeList.end() ) // not last merge
			{
				std::pair<int64_t,int64_t> nextMaster = next->masterBlocks;
				std::pair<int64_t,int64_t> nextSlave = next->slaveBlocks;

				if( curMaster.reversed ) this->reverseCoordinates( nextMaster, masterSize );
				if( curSlave.reversed ) this->reverseCoordinates( nextSlave, slaveSize );

				if( cur->masterId == next->masterId )
				{
					if( !curMaster.reversed ) fwdMerge = (cur->masterBlocks).first <= nextMaster.first;
					if( curMaster.reversed ) fwdMerge = (cur->masterBlocks).second <= nextMaster.second;
				}
				else // cur->slaveId == next->slaveId
				{
					if( !curSlave.reversed ) fwdMerge = (cur->slaveBlocks).first <= nextSlave.first;
					if( curSlave.reversed ) fwdMerge = (cur->slaveBlocks).second <= nextSlave.second;
				}

				// WARNING => da SISTEMARE!
				if( !curAlignFailed && !prevAlignFailed && fwdMerge != prevFwdMerge )
				{
					std::cerr << "nextSlave " << next->slaveId << ":\talign=(" << nextSlave.first << "," << nextSlave.second << ")"
					<< " reversed*=" << (next->align)->isCtgReversed() << " size=" << slaveSize << " hom=" << (next->align)->main_homology() << std::endl;

					std::cerr << "cambio direzione di merge da " << (fwdMerge ? "fwd" : "rev") << " a " << (prevFwdMerge ? "fwd" : "rev")
					<< ": master=" << curMaster.id << "=>" << prevMaster.id << " slave=" << curSlave.id << "=>" << prevSlave.id
					<< "\talign=" << (cur->align)->main_homology()
					<< " failed=" << this->isAlignmentFailed( *(cur->align), cur->masterPos, cur->slavePos, curMaster.reversed ) << std::endl;
					return;
				}
			}

			if( !curAlignFailed && !prevAlignFailed ) // good alignments (it's possible to fill the gap between two master contigs)
			{
				if( fwdMerge )
				{
					if( prevMaster.id == curMaster.id )
					{
						if( prevMaster.pos > (cur->masterMainAlign).first ) // current blocks-region overlapping previous one
						{
							if( prevMaster.pos > (cur->masterMainAlign).second ) // current blocks-region totally included in the previous one
							{
								curAlignFailed = true;
								curMaster.pos = prevMaster.pos;
								curSlave.pos = prevSlave.pos;
							}
							else
							{
								this->appendMasterToPctg( pctg, curMaster, prevMaster.pos, cur->masterMainAlign.second );
								curMaster.pos = cur->masterMainAlign.second + 1;
								curSlave.pos = cur->slaveMainAlign.second + 1;
							}
						}
						else // current blocks-region NOT overlapping previous one
						{
							// da dove ero arrivato fino a regione dei blocchi (compresa)
							this->appendMasterToPctg( pctg, curMaster, prevMaster.pos, cur->masterMainAlign.first - 1 );

							this->appendBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

							curMaster.pos = cur->masterMainAlign.second + 1;
							curSlave.pos = cur->slaveMainAlign.second + 1;
						}

					}
					else
					{
						if( prevSlave.pos > (cur->slaveMainAlign).first ) // current blocks-region overlapping previous one
						{
							if( prevSlave.pos > (cur->slaveMainAlign).second ) // current blocks-region totally included in the previous one
							{
								std::cerr << "totally included block found!" << std::endl;
								if( next != mergeList.end() && curMaster.id == next->masterId )
								{
									// completo pctg corrente con il la parte mancante del master precedente
									this->appendMasterToPctg( pctg, prevMaster, prevMaster.pos, prevMaster.limits.second );
									pctgList.push_back(pctg);

									// considero un nuovo pctg costituito dal master dall'inizio fino alla regione dei blocchi (compresa)
									pctg = PairedContig();

									this->appendMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterMainAlign.first-1 );

									this->appendBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
																	curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

									curMaster.pos = cur->masterMainAlign.second + 1;
									curSlave.pos = cur->slaveMainAlign.second + 1;
								}
								else
								{
									curMaster.pos = prevMaster.pos;
									curSlave.pos = prevSlave.pos;
								}

							}
							else // current blocks-region partially included in the previous one
							{
								// use slave to the end of the blocks-region
								this->appendSlaveToPctg( pctg, curSlave, prevSlave.pos, cur->slaveMainAlign.second);
								pctg.addMasterCtgId( curMaster.id );

								curMaster.pos = cur->masterMainAlign.second + 1;
								curSlave.pos = cur->slaveMainAlign.second + 1;
							}
						}
						else // current blocks-region NOT overlapping previous one
						{
							// regione che fa il salto
							this->appendSlaveToPctg( pctg, curSlave, prevSlave.pos, cur->slaveMainAlign.first-1);
							// regione dei blocchi

							this->appendBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );
							//this->appendMasterToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second );
							curMaster.pos = cur->masterMainAlign.second + 1;
							curSlave.pos = cur->slaveMainAlign.second + 1;
						}
					}
				}
				else // reveverse merge
				{
					if( prevMaster.id == curMaster.id )
					{
						if( prevMaster.pos < (cur->masterMainAlign).second )
						{
							if( prevMaster.pos < (cur->masterMainAlign).first ) // current blocks-region totally included in the previous one
							{
								curAlignFailed = true;
								curMaster.pos = prevMaster.pos;
								curSlave.pos = prevSlave.pos;
							}
							else
							{
								this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.first, prevMaster.pos );
								curMaster.pos = cur->masterMainAlign.first - 1;
								curSlave.pos = cur->slaveMainAlign.first - 1;
							}
						}
						else
						{
							// da dove ero arrivato fino a regione dei blocchi (compresa)
							this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.second+1, prevMaster.pos );
							this->prependBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															 curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

							curMaster.pos = cur->masterMainAlign.first - 1;
							curSlave.pos = cur->slaveMainAlign.first - 1;
						}
					}
					else
					{
						if( prevSlave.pos < (cur->slaveMainAlign).second )
						{
							if( prevSlave.pos < (cur->slaveMainAlign).first )
							{
								std::cerr << "totally included block found!" << std::endl;
								if( next != mergeList.end() && curMaster.id == next->masterId )
								{
									// dall'inizio fino a da dove ero arrivato
									this->prependMasterToPctg( pctg, prevMaster, prevMaster.limits.first, prevMaster.pos );
									pctgList.push_back(pctg);

									// dall'inizio fino alla regione dei blocchi (compresa)
									pctg = PairedContig();

									this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.second+1, cur->masterPos.second );

									this->prependBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
																	 curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

									curMaster.pos = cur->masterMainAlign.first - 1;
									curSlave.pos = cur->slaveMainAlign.first - 1;
								}
								else
								{
									curMaster.pos = prevMaster.pos;
									curSlave.pos = prevSlave.pos;
								}
							}
							else
							{
								// use slave to the end of the blocks-region
								this->prependSlaveToPctg( pctg, curSlave, cur->slaveMainAlign.first, prevSlave.pos );
								pctg.addMasterCtgId( curMaster.id );

								curMaster.pos = cur->masterMainAlign.first - 1;
								curSlave.pos = cur->slaveMainAlign.first - 1;
							}
						}
						else
						{
							// regione che fa il salto
							this->prependSlaveToPctg( pctg, curSlave, cur->slaveMainAlign.second + 1, prevSlave.pos );

							// regione dei blocchi
							this->prependBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															 curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );
							//this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second );
							curMaster.pos = cur->masterMainAlign.first - 1;
							curSlave.pos = cur->slaveMainAlign.first - 1;
						}
					}
				}
			}
			else // bad current alignment or bad previous alignment
			{
				// estendi il pctg
				if( fwdMerge )
				{
					if( prevMaster.id == curMaster.id )
					{
						if( !curAlignFailed ) // current alignment good, previous alignment failed
						{
							// da dove ero arrivato fino a regione dei blocchi (compresa)
							this->appendMasterToPctg( pctg, curMaster, prevMaster.pos, cur->masterMainAlign.first-1 );
							pctg.addMasterCtgId( curMaster.id );

							this->appendBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

							curMaster.pos = cur->masterMainAlign.second + 1;
							curSlave.pos = cur->slaveMainAlign.second + 1;
						}
						else
						{
							curMaster.pos = prevMaster.pos;
							curSlave.pos = prevSlave.pos;
						}
					}
					else // master contigs linked by slave
					{
						if( !curAlignFailed )
						{
							// completo pctg corrente con il la parte mancante del master precedente
							this->appendMasterToPctg( pctg, prevMaster, prevMaster.pos, prevMaster.limits.second );
							pctgList.push_back(pctg);

							// considero un nuovo pctg costituito dal master dall'inizio fino alla regione dei blocchi (compresa)
							pctg = PairedContig();

							this->appendMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterMainAlign.first-1 );

							this->appendBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

							curMaster.pos = cur->masterMainAlign.second + 1;
							curSlave.pos = cur->slaveMainAlign.second + 1;
						}
						else
						{
							// completo pctg corrente con il la parte mancante del master precedente
							this->appendMasterToPctg( pctg, prevMaster, prevMaster.pos, prevMaster.limits.second );
							pctg.addMasterCtgId( prevMaster.id );
							pctgList.push_back(pctg);

							// considero un nuovo pctg
							pctg = PairedContig();
							// dall'inizio fino alla regione dei blocchi (compresa)
							//this->appendMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterMainAlign.second );
							curMaster.pos = cur->masterPos.first;
							curSlave.pos = cur->slavePos.first;
						}
					}
				}
				else
				{
					if( prevMaster.id == curMaster.id )
					{
						if( !curAlignFailed )
						{
							// aggiungo master dalla regione dei blocchi (compresa) fino a dove ero arrivato
							this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.second+1, prevMaster.pos );
							pctg.addMasterCtgId( curMaster.id );

							this->prependBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															 curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

							curMaster.pos = cur->masterMainAlign.first - 1;
							curSlave.pos = cur->slaveMainAlign.first - 1;
						}
						else
						{
							curMaster.pos = prevMaster.pos;
							curSlave.pos = prevSlave.pos;
						}
					}
					else
					{
						if( !curAlignFailed )
						{
							// dall'inizio fino a da dove ero arrivato
							this->prependMasterToPctg( pctg, prevMaster, prevMaster.limits.first, prevMaster.pos );
							pctg.addMasterCtgId( prevMaster.id );
							pctgList.push_back(pctg);

							// dall'inizio fino alla regione dei blocchi (compresa)
							pctg = PairedContig();

							this->prependMasterToPctg( pctg, curMaster, cur->masterMainAlign.second+1, cur->masterPos.second );
							pctg.addMasterCtgId( curMaster.id );

							this->prependBlocksRegionToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterMainAlign.second,
															 curSlave, cur->slaveMainAlign.first, cur->slaveMainAlign.second );

							curMaster.pos = cur->masterMainAlign.first - 1;
							curSlave.pos = cur->slaveMainAlign.first - 1;
						}
						else
						{
							// dall'inizio fino a da dove ero arrivato
							this->prependMasterToPctg( pctg, prevMaster, prevMaster.limits.first, prevMaster.pos );
							pctg.addMasterCtgId( prevMaster.id );
							pctgList.push_back(pctg);

							pctg = PairedContig();
							// dalla regione dei blocchi (compresa) fino alla fine
							//this->appendMasterToPctg( pctg, curMaster, cur->masterMainAlign.first, cur->masterPos.second );
							curMaster.pos = cur->masterPos.second;
							curSlave.pos = cur->slavePos.second;
						}
					}
				}

				//prevAlignFailed = curAlignFailed;
			}
		}

		if( next == mergeList.end() ) // complete last merge
		{
			//bool curAlignFailed = this->isAlignmentFailed( *(cur->align), cur->masterPos, cur->slavePos, curMaster.reversed );

			if( !curAlignFailed && !prevAlignFailed ) // current alignment and previous alignment are good
			{
				if( fwdMerge )
				{
					if( (cur->masterPos).second - (cur->masterMainAlign).second >= (cur->slavePos).second - (cur->slaveMainAlign).second )
					{
						this->appendMasterToPctg( pctg, curMaster, curMaster.pos, cur->masterPos.second );
						pctg.addMasterCtgId( curMaster.id );
					}
					else
					{
						this->appendSlaveToPctg( pctg, curSlave, curSlave.pos, (cur->slavePos).second );
						pctg.addSlaveCtgId( curSlave.id );
					}
				}
				else
				{
					if( cur->masterMainAlign.first - cur->masterPos.first >= cur->slaveMainAlign.first - cur->slavePos.first )
					{
						this->prependMasterToPctg( pctg, curMaster, cur->masterPos.first, curMaster.pos );
						pctg.addMasterCtgId( curMaster.id );
					}
					else
					{
						this->prependSlaveToPctg( pctg, curSlave, cur->slavePos.first, curSlave.pos );
						pctg.addSlaveCtgId( curSlave.id );
					}
				}
			}
			else
			{
				if( fwdMerge )
				{
					if( prevMaster.id == curMaster.id )
					{
						if( !curAlignFailed )
						{
							//this->appendMasterToPctg( pctg, curMaster, curMaster.pos, (cur->masterMainAlign).second );
							//pctg.addMasterCtgId( curMaster.id );

							if( (cur->masterPos).second - (cur->masterMainAlign).second >= (cur->slavePos).second - (cur->slaveMainAlign).second )
							{
								this->appendMasterToPctg( pctg, curMaster, curMaster.pos, (cur->masterPos).second );
								pctg.addMasterCtgId( curMaster.id );
							}
							else
							{
								this->appendSlaveToPctg( pctg, curSlave, curSlave.pos, (cur->slavePos).second );
								pctg.addSlaveCtgId( curSlave.id );
							}
						}
						else
						{
							this->appendMasterToPctg( pctg, curMaster, curMaster.pos, (cur->masterPos).second );
							pctg.addMasterCtgId( curMaster.id );
						}
					}
					else
					{
						if( !curAlignFailed )
						{
							// dall'inizio fino alla regione dei blocchi (compresa)
							//this->appendMasterToPctg( pctg, curMaster, cur->masterPos.first, cur->masterMainAlign.second );

							if( (cur->masterPos).second - (cur->masterMainAlign).second >= (cur->slavePos).second - (cur->slaveMainAlign).second )
							{
								this->appendMasterToPctg( pctg, curMaster, curMaster.pos, (cur->masterPos).second );
								pctg.addMasterCtgId( curMaster.id );
							}
							else
							{
								this->appendSlaveToPctg( pctg, curSlave, curSlave.pos, (cur->slavePos).second );
								pctg.addSlaveCtgId( curSlave.id );
							}
						}
						else
						{
							this->appendMasterToPctg( pctg, curMaster, curMaster.pos, (cur->masterPos).second );
							pctg.addMasterCtgId( curMaster.id );
						}
					}
				}
				else
				{
					if( prevMaster.id == curMaster.id )
					{
						if( !curAlignFailed )
						{
							if( cur->masterMainAlign.first - cur->masterPos.first >= cur->slaveMainAlign.first - cur->slavePos.first)
							{
								this->prependMasterToPctg( pctg, curMaster, (cur->masterPos).first, curMaster.pos );
								pctg.addMasterCtgId( curMaster.id );
							}
							else
							{
								this->prependSlaveToPctg( pctg, curSlave, (cur->slavePos).first, curSlave.pos );
								pctg.addSlaveCtgId( curSlave.id );
							}
						}
						else
						{
							this->prependMasterToPctg( pctg, curMaster, (cur->masterPos).first, curMaster.pos );
							pctg.addMasterCtgId( curMaster.id );
						}
					}
					else
					{
						if( !curAlignFailed )
						{
							if( cur->masterMainAlign.first - cur->masterPos.first >= cur->slaveMainAlign.first - cur->slavePos.first)
							{
								this->prependMasterToPctg( pctg, curMaster, (cur->masterPos).first, curMaster.pos );
								pctg.addMasterCtgId( curMaster.id );
							}
							else
							{
								this->appendSlaveToPctg( pctg, curSlave, (cur->slavePos).first, curSlave.pos );
								pctg.addSlaveCtgId( curSlave.id );
							}
						}
						else
						{
							this->prependMasterToPctg( pctg, curMaster, (cur->masterPos).first, curMaster.pos );
							pctg.addMasterCtgId( curMaster.id );
						}
					}
				}
			}

			pctgList.push_back(pctg);
			pctg = PairedContig();
		}

		prevAlignFailed = curAlignFailed;

		if( prevMaster.id != -1 && prevMaster.id != curMaster.id ) delete prevMaster.ctg;
		if( prevSlave.id != -1 && prevSlave.id != curSlave.id ) delete prevSlave.ctg;

		prevMaster = curMaster;
		prevSlave = curSlave;

		prevFwdMerge = fwdMerge;

	} // end of while

	if( curMaster.id != -1 ) delete curMaster.ctg;
	if( curSlave.id != -1 ) delete curSlave.ctg;
}

bool PctgBuilder::isAlignmentFailed( BestCtgAlignment &align, std::pair<uint64_t,uint64_t> masterPos, std::pair<uint64_t,uint64_t> slavePos, bool isMasterRev )
{
	if( align.main_homology() < MIN_HOMOLOGY ) return true;

	size_t masterSize = align[0].a_size();
	size_t slaveSize = align[0].b_size();

	if( !isMasterRev )
	{
		if( align.left().homology() < MIN_HOMOLOGY && masterPos.first == 0 && slavePos.first == 0 ) return true;
		if( align.right().homology() < MIN_HOMOLOGY && masterPos.second == masterSize-1 && slavePos.second == slaveSize-1 ) return true;
	}
	else
	{
		if( align.right().homology() < MIN_HOMOLOGY && masterPos.first == 0 && slavePos.first == 0 ) return true;
		if( align.left().homology() < MIN_HOMOLOGY && masterPos.second == masterSize-1 && slavePos.second == slaveSize-1 ) return true;
	}

	return false;
}


void PctgBuilder::appendMasterToPctg( PairedContig &pctg, int32_t id, Contig &ctg, int32_t start, int32_t end, bool rev )
{
	if( end < start || start < 0 || end >= ctg.size() ) return;

	pctg.addMasterCtgId(id);

	int64_t idx = pctg.size();
	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg.at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_back( CtgInPctgInfo( id, start, end, rev, true ) );
}

void PctgBuilder::appendSlaveToPctg( PairedContig &pctg, int32_t id, Contig &ctg, int32_t start, int32_t end, bool rev )
{
	if( end < start || start < 0 || end >= ctg.size() ) return;

	pctg.addSlaveCtgId(id);

	int64_t idx = pctg.size();
	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg.at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_back( CtgInPctgInfo( id, start, end, rev, false ) );
}

void PctgBuilder::appendBlocksRegionToPctg( PairedContig &pctg, int32_t m_id, Contig &m_ctg, int32_t m_start, int32_t m_end, bool m_rev,
											int32_t s_id, Contig &s_ctg, int32_t s_start, int32_t s_end, bool s_rev )
{
	pctg.addMasterCtgId(m_id);
	pctg.addSlaveCtgId(s_id);

	int64_t master_int = (m_end >= m_start) ? m_end - m_start + 1 : 0;
	int64_t slave_int = (s_end >= s_start) ? s_end - s_start + 1 : 0;

	int64_t large_int = (master_int >= slave_int) ? master_int : slave_int;
	int64_t small_int = (master_int >= slave_int) ? slave_int : master_int;

	// if regions are similar in length append master's blocks region and return
	if( small_int >= 0.97*large_int ) return appendMasterToPctg( pctg, m_id, m_ctg, m_start, m_end, m_rev );

	std::vector<double> masterScore, slaveScore;

	masterScore = computeZScore( masterBam, m_id, m_start, m_end );
	slaveScore = computeZScore( slaveBam, s_id, s_start, s_end );

	size_t masterEvid = 0;
	size_t slaveEvid = 0;

	for( size_t i=0; i < masterScore.size(); i++ )
	{
		double m_score = masterScore[i]; if( m_score < 0 ) m_score = -m_score;
		double s_score =  slaveScore[i]; if( s_score < 0 ) s_score = -s_score;

		if( s_score < m_score && s_score != 0 ) slaveEvid++; else if( s_score < m_score ) masterEvid++;
		if( m_score < s_score && m_score != 0 ) masterEvid++; else if( m_score < s_score ) slaveEvid++;
	}

	if( masterEvid >= slaveEvid ) return appendMasterToPctg( pctg, m_id, m_ctg, m_start, m_end, m_rev );

	return appendSlaveToPctg( pctg, s_id, s_ctg, s_start, s_end, s_rev );
}


void PctgBuilder::buildPctgs( std::list<PairedContig> &pctgList, MergeBlockLists &mergeLists )
{
	for( MergeBlockLists::iterator it = mergeLists.begin(); it != mergeLists.end(); ++it )
	{
		if( it->size() == 0 ) continue;
		this->buildPctgs( pctgList, *it );
	}
}


void PctgBuilder::buildPctgs( std::list<PairedContig> &pctgList, std::list<MergeBlock> &ml )
{
	if( ml.size() == 0 )
	{
		std::cerr << "buildPctgs: MergeBlocks' list empty. This shouldn't happen!" << std::endl;
		return;
	}

	std::list<MergeBlock>::iterator it, mb, mb_next;

	PairedContig pctg;

	int32_t m_pos = 0;
	int32_t s_pos = 0;

	Contig *master_ctg, *slave_ctg;
	int32_t prev_mid, prev_sid;

	it = ml.begin();
	while( it != ml.end() )
	{
		mb = it;
		mb_next = ++it;

		if( mb == ml.begin() ) // first merge
		{
			master_ctg = new Contig( this->loadMasterContig( mb->m_id ) );
			slave_ctg = new Contig( this->loadSlaveContig( mb->s_id ) );

			if( mb->m_rev ) reverse_complement(*master_ctg);
			if( mb->s_rev ) reverse_complement(*slave_ctg);

			// add first tail
			int32_t m_tail = ( mb->m_ltail ) ? mb->m_start : 0;
			int32_t s_tail = 0; //( mb->ext_slave_prev && mb->s_ltail ) ? mb->s_start : 0;

			if( m_tail >= s_tail && m_tail > 0 ) this->appendMasterToPctg( pctg, mb->m_id, *master_ctg, 0, mb->m_start - 1, mb->m_rev );
			if( s_tail > m_tail && s_tail > 0 ) this->appendSlaveToPctg( pctg, mb->s_id, *slave_ctg, 0, mb->s_start - 1, mb->s_rev );

			// add block region
			this->appendBlocksRegionToPctg( pctg, mb->m_id, *master_ctg, mb->m_start, mb->m_end, mb->m_rev, mb->s_id, *slave_ctg, mb->s_start, mb->s_end, mb->s_rev );
		}
		else // not first block
		{
			if( mb->m_id == prev_mid )
			{
				delete slave_ctg;
				slave_ctg = new Contig( this->loadSlaveContig( mb->s_id ) );
				if( mb->s_rev ) reverse_complement(*slave_ctg);

				if( m_pos <= mb->m_start )
				{
					// fill the gap
					this->appendMasterToPctg( pctg, mb->m_id, *master_ctg, m_pos, mb->m_start - 1, mb->m_rev);
					// add block region
					this->appendBlocksRegionToPctg( pctg, mb->m_id, *master_ctg, mb->m_start, mb->m_end, mb->m_rev, mb->s_id, *slave_ctg, mb->s_start, mb->s_end, mb->s_rev );
				}
				else // current merge block overlaps previous one
				{
					this->appendMasterToPctg( pctg, mb->m_id, *master_ctg, m_pos, mb->m_end, mb->m_rev);
				}
			}
			else // mb->s_id == prev_sid
			{
				delete master_ctg;
				master_ctg = new Contig( this->loadMasterContig( mb->m_id ) );
				if( mb->m_rev ) reverse_complement(*master_ctg);

				if( s_pos <= mb->s_start )
				{
					// fill the gap
					this->appendSlaveToPctg( pctg, mb->s_id, *slave_ctg, s_pos, mb->s_start - 1, mb->s_rev);
					// add block region
					this->appendBlocksRegionToPctg( pctg, mb->m_id, *master_ctg, mb->m_start, mb->m_end, mb->m_rev, mb->s_id, *slave_ctg, mb->s_start, mb->s_end, mb->s_rev );
				}
				else
				{
					this->appendSlaveToPctg( pctg, mb->s_id, *slave_ctg, s_pos, mb->s_end, mb->s_rev );
					pctg.addMasterCtgId( mb->m_id );
				}
			}
		}

		if( mb_next == ml.end() ) // for last block, add tail if possible
		{
			int32_t m_size = master_ctg->size();
			int32_t s_size = slave_ctg->size();

			int32_t m_tail = ( mb->m_rtail ) ? m_size - mb->m_end - 1 : 0;
			int32_t s_tail = 0; //( mb->ext_slave_next && mb->s_rtail ) ? s_size - mb->s_end - 1 : 0;

			if( m_tail >= s_tail && m_tail > 0 ) this->appendMasterToPctg( pctg, mb->m_id, *master_ctg, mb->m_end+1, m_size-1, mb->m_rev );
			if( s_tail > m_tail && s_tail > 0 ) this->appendSlaveToPctg( pctg, mb->s_id, *slave_ctg, mb->s_end+1, s_size-1, mb->s_rev );

			delete master_ctg;
			delete slave_ctg;
		}

		prev_mid = mb->m_id;
		prev_sid = mb->s_id;

		m_pos = mb->m_end+1;
		s_pos = mb->s_end+1;
	}

	if( pctg.size() > 0 ) pctgList.push_back(pctg);
}


void PctgBuilder::splitMergeBlocksByInclusions( MergeBlockLists &ml_in )
{
	const RefVector& master_ref = masterBam.GetReferenceData();
	const RefVector& slave_ref = slaveBam.GetReferenceData();

	MergeBlockLists tmp;
	MergeBlock mb_cur, mb_next, mb_prev;

	for( MergeBlockLists::iterator ml = ml_in.begin(); ml != ml_in.end(); ++ml )
	{
		bool first = true;
		bool split_prev = false;
		bool fwd_merge = true;
		bool fwd_merge_prev = true;

		bool master_rev, slave_rev;

		std::list<MergeBlock>::iterator mb_it;

		std::list<MergeBlock> ml_new;

		mb_it = ml->begin();
		while( mb_it != ml->end() )
		{
			mb_cur = *(mb_it++);
			if( mb_it != ml->end() ) mb_next = *(mb_it);

			if( first ) // first
			{
				first = false;
				//mb_cur.m_rev = false;
				//mb_cur.s_rev = mb_cur.align_rev;

				// update mergeblock coordinates wrt master ctg orientation
				if( mb_cur.m_rev )
				{
					int32_t m_size = master_ref.at( mb_cur.m_id ).RefLength;
					int32_t tmp_start = mb_cur.m_start;
					mb_cur.m_start = m_size - mb_cur.m_end - 1;
					mb_cur.m_end = m_size - tmp_start - 1;

					std::swap( mb_cur.m_ltail, mb_cur.m_rtail );
				}

				// update mergeblock coordinates wrt slave ctg orientation
				if( mb_cur.s_rev )
				{
					int32_t s_size = slave_ref.at( mb_cur.s_id ).RefLength;
					int32_t tmp_start = mb_cur.s_start;
					mb_cur.s_start = s_size - mb_cur.s_end - 1;
					mb_cur.s_end = s_size - tmp_start - 1;

					std::swap( mb_cur.s_ltail, mb_cur.s_rtail );
				}

				first = false;
				ml_new.push_back(mb_cur);

				mb_prev = mb_cur;
			}
			else // not first
			{
				// update mergeblock coordinates wrt master ctg orientation
				if( mb_cur.m_rev )
				{
					int32_t m_size = master_ref.at( mb_cur.m_id ).RefLength;
					int32_t tmp_start = mb_cur.m_start;
					mb_cur.m_start = m_size - mb_cur.m_end - 1;
					mb_cur.m_end = m_size - tmp_start - 1;

					std::swap( mb_cur.m_ltail, mb_cur.m_rtail );
				}

				// update mergeblock coordinates wrt slave ctg orientation
				if( mb_cur.s_rev )
				{
					int32_t s_size = slave_ref.at( mb_cur.s_id ).RefLength;
					int32_t tmp_start = mb_cur.s_start;
					mb_cur.s_start = s_size - mb_cur.s_end - 1;
					mb_cur.s_end = s_size - tmp_start - 1;

					std::swap( mb_cur.s_ltail, mb_cur.s_rtail );
				}

				if( mb_prev.m_id == mb_cur.m_id ) // jumped from master
				{
					// previous merge-block fully included in current one
					if( mb_prev.m_start > mb_cur.m_start && mb_prev.m_end <= mb_cur.m_end )
					{
						//std::cerr << "INCLUDED blocks in master ctg id = " << mb_cur.m_id << std::endl;
						while( ml_new.size() > 0 && ml_new.back().m_start > mb_cur.m_start && ml_new.back().m_end <= mb_cur.m_end && ml_new.back().m_id == mb_cur.m_id )
							ml_new.pop_back();

						if( ml_new.size() > 0 && ml_new.back().m_id != mb_cur.m_id && ml_new.back().s_id != mb_cur.s_id )
						{
							ml_new.back().ext_slave_next = false;
							tmp.push_back( ml_new );

							ml_new.clear();
						}

						ml_new.push_back(mb_cur);
						mb_prev = mb_cur;
					}
					else if( mb_prev.m_start > mb_cur.m_start ) // previous merge-block begins after current one
					{
						/*if( mb_it != ml->end() && mb_next.s_id == mb_cur.s_id )
						{
							ml_new.back().ext_slave_next = false;
							tmp.push_back( ml_new );

							ml_new.clear();
							mb_it->ext_slave_prev = false;

							continue;
						}*/

						ml_new.back().ext_slave_next = false;
						break;
					}
					else if( mb_prev.m_end >= mb_cur.m_end ) // current merge-block fully included in previous one
					{
						if( mb_it != ml->end() ) // if exists a next MergeBlock
						{
							// jumping to next with master
							if( mb_cur.m_id == mb_next.m_id ) continue;

							// jumping to next with slave
							ml_new.back().ext_slave_next = false;
							tmp.push_back( ml_new );

							ml_new.clear();
							mb_it->ext_slave_prev = false; // next

							first = true;
						}

						continue;
					}
					else // current merge-block not included in the previous one
					{
						ml_new.push_back(mb_cur);
						mb_prev = mb_cur;
					}
				}
				else // jumped from slave
				{
					// previous merge-block fully included in current one
					if( mb_prev.s_start > mb_cur.s_start && mb_prev.s_end <= mb_cur.s_end )
					{
						//std::cerr << "INCLUDED blocks in master ctg id = " << mb_cur.m_id << std::endl;
						while( ml_new.size() > 0 && ml_new.back().s_start > mb_cur.s_start && ml_new.back().s_end <= mb_cur.s_end && ml_new.back().s_id == mb_cur.s_id )
							ml_new.pop_back();

						if( ml_new.size() > 0 && ml_new.back().m_id != mb_cur.m_id && ml_new.back().s_id != mb_cur.s_id )
						{
							ml_new.back().ext_slave_next = false;
							tmp.push_back( ml_new );

							ml_new.clear();
						}

						ml_new.push_back(mb_cur);
						mb_prev = mb_cur;
					}
					else if( mb_prev.s_start > mb_cur.s_start ) // previous merge-block begins after current one
					{
						/*if( mb_it != ml->end() && mb_next.m_id == mb_cur.m_id )
						{
							ml_new.back().ext_slave_next = false;
							tmp.push_back( ml_new );

							ml_new.clear();
							mb_it->ext_slave_prev = false;

							continue;
						}*/

						ml_new.back().ext_slave_next = false;
						break;
					}
					else if( mb_prev.s_end >= mb_cur.s_end ) // current merge-block fully included in previous one
					{
						if( mb_it != ml->end() ) // if exists a next MergeBlock
						{
							// jumping to next with slave
							if( mb_cur.s_id == mb_next.s_id ) continue;

							// jumping to next with master
							ml_new.back().ext_slave_next = false;
							tmp.push_back( ml_new );

							ml_new.clear();
							mb_next.ext_slave_prev = false;

							first = true;
						}

						continue;
					}
					else // current merge-block not included in the previous one
					{
						ml_new.push_back(mb_cur);
						mb_prev = mb_cur;
					}
				}
			}
		} // end of while (list processing)

		if( ml_new.size() > 0 ) tmp.push_back(ml_new);
	} // end of for (list of mergeblocks lists processing)

	ml_in = tmp;
}


void PctgBuilder::sortMergeBlocksByDirection( MergeBlockLists &ml )
{
	for( MergeBlockLists::iterator ml_in = ml.begin(); ml_in != ml.end(); ml_in++ )
	{
		std::list<MergeBlock>::iterator mb, first, second;

		if( ml_in->size() < 2 ) continue;

		mb = ml_in->begin();

		first = mb++;
		second = mb++;

		bool fwd_merge, master_rev, slave_rev;

		slave_rev = first->align_rev;

		if( first->m_id == second->m_id )
		{
			fwd_merge = first->m_start <= second->m_start;
		}
		else
		{
			if( !slave_rev ) fwd_merge = (first->s_start <= second->s_start);
			else fwd_merge = (first->s_start >= second->s_start);
		}

		if( !fwd_merge ) // reverse list
		{
			for( mb = ml_in->begin(); mb != ml_in->end(); ++mb ) std::swap( mb->ext_slave_next, mb->ext_slave_prev );
			ml_in->reverse();
		}
	}
}


void PctgBuilder::splitMergeBlocksByDirection( MergeBlockLists &ml_in )
{
	MergeBlockLists ml_out;
	int32_t master_id, slave_id;

	for( MergeBlockLists::iterator ml = ml_in.begin(); ml != ml_in.end(); ++ml )
	{
		bool first = true;
		bool split_prev = false;
		bool fwd_merge = true;
		bool fwd_merge_prev = true;

		bool master_rev, slave_rev;

		std::list<MergeBlock>::iterator mb, cur, next;

		std::list<MergeBlock> ml_new;

		mb = ml->begin();
		while( mb != ml->end() )
		{
			cur = mb;
			next = ++mb;

			if( first ) // first
			{
				master_id = cur->m_id;
				slave_id = cur->s_id;
				master_rev = false;
				slave_rev = cur->align_rev;

				cur->m_rev = master_rev;
				cur->s_rev = slave_rev;

				if( split_prev )
				{
					cur->ext_slave_prev = false;
					split_prev = false;
				}

				if( next != ml->end() ) // not last merge block
				{
					if( cur->m_id == next->m_id )
					{
						fwd_merge = cur->m_start <= next->m_start;
					}
					else // cur->s_id == next->s_id
					{
						if( !slave_rev ) fwd_merge = (cur->s_start <= next->s_start);
							else fwd_merge = (cur->s_start >= next->s_start);
					}
				}

				first = false;
				fwd_merge_prev = fwd_merge;
				ml_new.push_back(*cur);
			}
			else // not first
			{
				if( master_id == cur->m_id ) slave_rev = (master_rev && !cur->align_rev) || (!master_rev && cur->align_rev);
				if( slave_id == cur->s_id ) master_rev = (slave_rev && !cur->align_rev) || (!slave_rev && cur->align_rev);

				cur->m_rev = master_rev;
				cur->s_rev = slave_rev;

				if( next != ml->end() ) // not last merge block
				{
					if( cur->m_id == next->m_id )
					{
						if( !master_rev ) fwd_merge = (cur->m_start <= next->m_start);
							else fwd_merge = (cur->m_start >= next->m_start);
					}
					else // cur->s_id == next->s_id
					{
						if( !slave_rev ) fwd_merge = (cur->s_start <= next->s_start);
						else fwd_merge = (cur->s_start >= next->s_start);
					}

					if( fwd_merge != fwd_merge_prev )
					{
						if( ml_new.back().m_id == cur->m_id && cur->m_id == next->m_id )
						{
							//std::cerr << "merge changed \"direction\" when it shouldn't (3 successive blocks with same master_id)" << std::endl;

							ml_new.push_back(*cur);
							master_id = cur->m_id;
							slave_id = cur->s_id;
							continue;
						}

						if( ml_new.back().s_id == cur->s_id && cur->s_id == next->s_id )
						{
							//std::cerr << "merge changed \"direction\" when it shouldn't (3 successive blocks with same slave_id)" << std::endl;

							ml_new.push_back(*cur);
							master_id = cur->m_id;
							slave_id = cur->s_id;
							continue;
						}

						ml_new.back().ext_slave_next = false;
						split_prev = true;
						first = true;

						if( ml_new.size() > 0 ) ml_out.push_back(ml_new);
						ml_new.clear();
						continue;
					}
				}

				ml_new.push_back(*cur);

				master_id = cur->m_id;
				slave_id = cur->s_id;
			}
		} // end of while

		if( ml_new.size() > 0 ) ml_out.push_back(ml_new);
	}

	ml_in = ml_out;
}


void PctgBuilder::splitMergeBlocksByAlign( MergeBlockLists &ml_in )
{
	MergeBlockLists ml_out;

	for( MergeBlockLists::iterator ml = ml_in.begin(); ml != ml_in.end(); ++ml )
	{
		std::list<MergeBlock>::iterator mb, cur, next;

		std::list<MergeBlock> ml_new;

		mb = ml->begin();
		bool prev_failed = false;
		while( mb != ml->end() )
		{
			cur = mb;
			next = ++mb;

			if( !cur->align_ok ) // current merge block's alignment failed
			{
				prev_failed = true;
				continue;
			}

			if( prev_failed ) cur->ext_slave_prev = false;
			if( next != ml->end() && !next->align_ok ) cur->ext_slave_next = false;

			if( ml_new.size() > 0 )
			{
				if( !prev_failed && (ml_new.back().m_id == cur->m_id || ml_new.back().s_id == cur->s_id) )
				{
					ml_new.push_back(*cur);
				}
				else if( prev_failed && ml_new.back().m_id == cur->m_id )
				{
					ml_new.push_back(*cur);
				}
				else
				{
					ml_out.push_back(ml_new);
					ml_new.clear();
					ml_new.push_back(*cur);
				}
			}
			else
			{
				ml_new.push_back(*cur);
			}

			prev_failed = false;
		}

		// add last "merge list"
		if( ml_new.size() > 0 ) ml_out.push_back(ml_new);
	}

	ml_in = ml_out;
}


void PctgBuilder::alignMergeBlock( const CompactAssemblyGraph &graph, MergeBlock &mb ) const
{
	typedef CompactAssemblyGraph::Vertex Vertex;

	Vertex v = mb.vertex;
	const std::list<Block> &blocks_list = graph.getBlocks(v);

	const Block &firstBlock = blocks_list.front();
	const Block &lastBlock = blocks_list.back();

	const Frame &firstMasterFrame = firstBlock.getMasterFrame();
	const Frame &firstSlaveFrame = firstBlock.getSlaveFrame();
	const Frame &lastMasterFrame = lastBlock.getMasterFrame();
	const Frame &lastSlaveFrame = lastBlock.getSlaveFrame();

	int32_t masterStart = std::min( firstMasterFrame.getBegin(), lastMasterFrame.getBegin() );
	int32_t masterEnd = std::max( firstMasterFrame.getEnd(), lastMasterFrame.getEnd() );
	int32_t slaveStart = std::min( firstSlaveFrame.getBegin(), lastSlaveFrame.getBegin() );
	int32_t slaveEnd = std::max( firstSlaveFrame.getEnd(), lastSlaveFrame.getEnd() );

	// load contig that should be merged with pctg.
	Contig masterCtg = this->loadMasterContig(mb.m_id);
	Contig slaveCtg = this->loadSlaveContig(mb.s_id);

	// find best alignment between the contigs
	BestCtgAlignment *bestAlign = new BestCtgAlignment();
	this->findBestAlignment( *bestAlign, masterCtg, masterStart, masterEnd, slaveCtg, slaveStart, slaveEnd, blocks_list );

	// find start/end positions of the best alignment
	std::pair<uint64_t,uint64_t> alignStart, alignEnd, alignStartTmp, alignEndTmp;

	mb.align_ok = true;

	if( bestAlign->main_homology() >= MIN_HOMOLOGY )
	{
		first_match_pos( bestAlign->at(0), alignStart );
		last_match_pos( bestAlign->at(bestAlign->size()-1), alignEnd );

		// compute masterCtg tails (i1,i2) and slaveCtg tails (j1,j2)
		uint64_t i1 = alignStart.first;
		uint64_t i2 = masterCtg.size() - alignEnd.first - 1;
		uint64_t j1 = alignStart.second;
		uint64_t j2 = slaveCtg.size() - alignEnd.second - 1;

		const MyAlignment& left = bestAlign->left();
		const MyAlignment& right = bestAlign->right();

		size_t mt = 0.3 * masterCtg.size();
		size_t st = 0.3 * slaveCtg.size();

		uint64_t left_min_len = 0.7 * std::min(i1,j1);
		uint64_t right_min_len = 0.7 * std::min(i2,j2);
		uint64_t threshold = std::min( size_t(100), std::min(mt,st) );

		bool s_ltail = bestAlign->isCtgReversed() ? mb.s_rtail : mb.s_ltail;
		bool s_rtail = bestAlign->isCtgReversed() ? mb.s_ltail : mb.s_rtail;

		if( mb.m_ltail && s_ltail && std::min(i1,j1) >= threshold )
		{
			if( this->is_good(left,left_min_len) )
			{
				first_match_pos(left,alignStart);
				if( bestAlign->is_left_rev() ) std::swap( alignStart.first, alignStart.second );
			}
			else
			{
				mb.align_ok = false;
			}
		}

		if( mb.m_rtail && s_rtail && std::min(i2,j2) >= threshold )
		{
			if( this->is_good(right,right_min_len) )
			{
				//alignEndTmp = alignEnd;
				last_match_pos(right,alignEndTmp);

				if( bestAlign->is_right_rev() )
				{
					std::swap( alignEndTmp.first, alignEndTmp.second );

					alignEnd.first = alignEndTmp.first;
					alignEnd.second += alignEndTmp.second+1;
				}
				else
				{
					alignEnd.first += alignEndTmp.first+1;
					alignEnd.second = alignEndTmp.second;
				}

				//alignEnd.first += alignEndTmp.first+1;
				//alignEnd.second += alignEndTmp.second+1;
			}
			else
			{
				mb.align_ok = false;
			}
		}
	}
	else // bad alignment between blocks
	{
		mb.align_ok = false;		
		return;
	}

	if( bestAlign->isCtgReversed() )
	{
		uint64_t tempPos = alignStart.second;
		alignStart.second = slaveCtg.size() - alignEnd.second - 1;
		alignEnd.second = slaveCtg.size() - tempPos - 1;
	}

	mb.align_rev = bestAlign->isCtgReversed();

	mb.m_start = alignStart.first;
	mb.m_end = alignEnd.first;
	mb.s_start = alignStart.second;
	mb.s_end = alignEnd.second;
}


void PctgBuilder::getNextContigs( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::set<int32_t> &masterCtgs, std::set<int32_t> &slaveCtgs ) const
{
    typedef CompactAssemblyGraph::Vertex Vertex;
    typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;
    
    OutEdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::out_edges(v,graph);
    
    int32_t mid = graph.getBlocks(v).front().getMasterId();
    int32_t sid = graph.getBlocks(v).front().getSlaveId();
    
    masterCtgs.insert(mid);
    slaveCtgs.insert(sid);
        
    for( OutEdgeIterator e = ebegin; e != eend; e++ )
    {
        Vertex v_next = boost::target(*e,graph);
        this->getNextContigs( graph, v_next, masterCtgs, slaveCtgs );
    }
}

void PctgBuilder::getPrevContigs( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::set<int32_t> &masterCtgs, std::set<int32_t> &slaveCtgs ) const
{
    typedef CompactAssemblyGraph::Vertex Vertex;
    typedef CompactAssemblyGraph::InEdgeIterator InEdgeIterator;
    
    InEdgeIterator ebegin,eend;
    boost::tie(ebegin,eend) = boost::in_edges(v,graph);
    
    int32_t mid = graph.getBlocks(v).front().getMasterId();
    int32_t sid = graph.getBlocks(v).front().getSlaveId();
    
    masterCtgs.insert(mid);
    slaveCtgs.insert(sid);
        
    for( InEdgeIterator e = ebegin; e != eend; e++ )
    {
        Vertex v_next = boost::source(*e,graph);
        this->getPrevContigs( graph, v_next, masterCtgs, slaveCtgs );
    }
}


bool PctgBuilder::solveForks( CompactAssemblyGraph &graph, std::vector<MergeBlock> &mbv ) const
{
	typedef CompactAssemblyGraph::Vertex Vertex;
	typedef CompactAssemblyGraph::VertexIterator VertexIterator;
	typedef CompactAssemblyGraph::Edge Edge;
	typedef CompactAssemblyGraph::EdgeIterator EdgeIterator;
	typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;
	typedef CompactAssemblyGraph::InEdgeIterator InEdgeIterator;
	
	VertexIterator vbegin,vend;
	EdgeIterator ebegin,eend;
	
	// merge-blocks initialization
	boost::tie(vbegin,vend) = boost::vertices(graph);
	for( VertexIterator v = vbegin; v != vend; v++ )
	{
		const std::list<Block> &blocks_list = graph.getBlocks(*v);
		
		mbv[*v].valid = true;
		
		mbv[*v].vertex = *v;
		mbv[*v].m_id = blocks_list.front().getMasterId();
		mbv[*v].s_id = blocks_list.front().getSlaveId();
				
		mbv[*v].ext_slave_next = true;
		mbv[*v].ext_slave_prev = true;
		
		mbv[*v].m_ltail = true;
		mbv[*v].m_rtail = true;
		mbv[*v].s_ltail = true;
		mbv[*v].s_rtail = true;
	}
	
	// handle putative repeats (nodes with io degree >= 2)
	boost::tie(vbegin,vend) = boost::vertices(graph);
	for( VertexIterator v = vbegin; v != vend; v++ ) 
	{
		int32_t in_deg = boost::in_degree(*v,graph);
		int32_t out_deg = boost::out_degree(*v,graph);
		
		if( in_deg >= 2 && out_deg >= 2 )
		{
			mbv[*v].valid = false;
			
			Vertex mv1, mv2, sv1, sv2;
			EdgeProperty edge_prop;
			
			double mw = 1.0;
			double sw = 1.0;
			
			InEdgeIterator ie_begin,ie_end;
			boost::tie(ie_begin,ie_end) = boost::in_edges(*v,graph);
			
			while( ie_begin != ie_end ) // remove input edges
			{
				edge_prop = boost::get(boost::edge_kind_t(), graph, *ie_begin);
				
				if( edge_prop.kind == MASTER_EDGE )
				{
					mv1 = boost::source(*ie_begin,graph); 
					mw = std::min(edge_prop.weight,mw);
				}
				else
				{
					sv1 = boost::source(*ie_begin,graph);
					sw = std::min(edge_prop.weight,sw);
				}
				
				boost::remove_edge( boost::source(*ie_begin,graph), *v, graph );
				boost::tie(ie_begin,ie_end) = boost::in_edges(*v,graph);
			}
			
			OutEdgeIterator oe_begin,oe_end;
			boost::tie(oe_begin,oe_end) = boost::out_edges(*v,graph);
			
			while( oe_begin != oe_end ) // remove output edges
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *oe_begin);
				if( edge_prop.kind == MASTER_EDGE )
				{
					mv2 = boost::target(*oe_begin,graph); 
					mw = std::min(edge_prop.weight,mw);
				}
				else 
				{
					sv2 = boost::target(*oe_begin,graph);
					sw = std::min(edge_prop.weight,sw);
				}
				
				boost::remove_edge( *v, boost::target(*oe_begin,graph), graph );
				boost::tie(oe_begin,oe_end) = boost::out_edges(*v,graph);
			}
			
			Edge e;
			
			edge_prop.kind = MASTER_EDGE; edge_prop.weight = mw;
			e = boost::add_edge( mv1, mv2, graph ).first;
			boost::put( boost::edge_kind_t(), graph, e, edge_prop );
			
			edge_prop.kind = SLAVE_EDGE; edge_prop.weight = sw;
			e = boost::add_edge( sv1, sv2, graph ).first;
			boost::put( boost::edge_kind_t(), graph, e, edge_prop );
		}
	}
	
	// handle bifurcations (putative mis-assemblies)
	boost::tie(vbegin,vend) = boost::vertices(graph);
	for( VertexIterator v = vbegin; v != vend; v++ )
	{
		const std::list<Block> &blocks_list = graph.getBlocks(*v);
		int32_t masterStart = std::min( blocks_list.front().getMasterFrame().getBegin(), blocks_list.back().getMasterFrame().getBegin() );
		int32_t slaveStart = std::min( blocks_list.front().getSlaveFrame().getBegin(), blocks_list.back().getSlaveFrame().getBegin() );
		
		int32_t in_deg = boost::in_degree(*v,graph);
		int32_t out_deg = boost::out_degree(*v,graph);
		
		if( in_deg < 2 && out_deg < 2 ) continue;
		
		Vertex mv, sv, ov; // master/slave linked vertices
		
		if( in_deg == 2 )
		{
			InEdgeIterator ebegin,eend, me, se;
			boost::tie(ebegin,eend) = boost::in_edges(*v,graph);
			
			double mw,sw;

			for( InEdgeIterator e = ebegin; e != eend; e++ )
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e);
				
				if( edge_prop.kind == MASTER_EDGE )
				{
					me = e;
					mv = boost::source(*e,graph); 
					mw = edge_prop.weight;
				}
				else
				{
					se = e;
					sv = boost::source(*e,graph);
					sw = edge_prop.weight;
				}
			}
			
			if( mw < 0 || sw < 0 ) continue;
			
			double w_diff = mw >= sw ? mw-sw : sw-mw;
			ForkType fork_type = w_diff >= 0.8 ? (mw >= sw ? MIS_SLAVE : MIS_MASTER) : UNKNOWN;
			
			if( fork_type == UNKNOWN ) continue;
			
			const std::list<Block> &masterNextBlocks = graph.getBlocks(mv);
			const std::list<Block> &slaveNextBlocks = graph.getBlocks(sv);

			int32_t nextMasterStart = std::min( masterNextBlocks.front().getMasterFrame().getBegin(), masterNextBlocks.back().getMasterFrame().getBegin() );
			int32_t nextSlaveStart = std::min( slaveNextBlocks.front().getSlaveFrame().getBegin(), slaveNextBlocks.back().getSlaveFrame().getBegin() );

			bool sharedFirstInMaster = masterStart <= nextMasterStart;
			bool sharedFirstInSlave = slaveStart <= nextSlaveStart;

			if( fork_type == MIS_MASTER )
			{
				std::cerr << "[debug] Found MASTER mis-assembly in ctg "  << blocks_list.front().getMasterId()
					<< " mw=" << mw << " sw=" << sw << " w_diff=" << w_diff << std::endl;

				if( sharedFirstInMaster ){ mbv[*v].m_rtail = false; mbv[mv].m_ltail = false; }
				else { mbv[*v].m_ltail = false; mbv[mv].m_rtail = false; }

				boost::remove_edge( *me, graph );
			}
			else if( fork_type == MIS_SLAVE )
			{
				std::cerr << "[debug] Found SLAVE mis-assembly in ctg "  << blocks_list.front().getSlaveId()
					<< " mw=" << mw << " sw=" << sw << " w_diff=" << w_diff << std::endl;

				if( sharedFirstInSlave ){ mbv[*v].s_rtail = false; mbv[sv].s_ltail = false; }
				else { mbv[*v].s_ltail = false; mbv[sv].s_rtail = false; }

				boost::remove_edge( *se, graph );
			}
		}
		
		if( out_deg == 2 )
		{
			OutEdgeIterator ebegin,eend, me, se;
			boost::tie(ebegin,eend) = boost::out_edges(*v,graph);
			
			double mw,sw;

			for( OutEdgeIterator e = ebegin; e != eend; e++ )
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e);
				
				if( edge_prop.kind == MASTER_EDGE )
				{
					me = e;
					mv = boost::target(*e,graph);
					mw = edge_prop.weight;
				}
				else
				{
					se = e;
					sv = boost::target(*e,graph);
					sw = edge_prop.weight;
				}
			}
			
			if( mw < 0 || sw < 0 ) continue;
			
			double w_diff = mw >= sw ? mw-sw : sw-mw;
			ForkType fork_type = w_diff >= 0.8 ? (mw >= sw ? MIS_SLAVE : MIS_MASTER) : UNKNOWN;
			
			if( fork_type == UNKNOWN ) continue;
			
			const std::list<Block> &masterNextBlocks = graph.getBlocks(mv);
			const std::list<Block> &slaveNextBlocks = graph.getBlocks(sv);

			int32_t nextMasterStart = std::min( masterNextBlocks.front().getMasterFrame().getBegin(), masterNextBlocks.back().getMasterFrame().getBegin() );
			int32_t nextSlaveStart = std::min( slaveNextBlocks.front().getSlaveFrame().getBegin(), slaveNextBlocks.back().getSlaveFrame().getBegin() );

			bool sharedFirstInMaster = masterStart <= nextMasterStart;
			bool sharedFirstInSlave = slaveStart <= nextSlaveStart;

			if( fork_type == MIS_MASTER )
			{
				std::cerr << "[debug] Found MASTER misassembly in ctg "  << blocks_list.front().getMasterId()
					<< " mw=" << mw << " sw=" << sw << " w_diff=" << w_diff << std::endl;

				if( sharedFirstInMaster ){ mbv[*v].m_rtail = false; mbv[mv].m_ltail = false; }
				else { mbv[*v].m_ltail = false; mbv[mv].m_rtail = false; }

				boost::remove_edge( me, graph );
			}
			else if( fork_type == MIS_SLAVE )
			{
				std::cerr << "[debug] Found SLAVE misassembly in ctg "  << blocks_list.front().getSlaveId()
					<< " mw=" << mw << " sw=" << sw << " w_diff=" << w_diff << std::endl;

				if( sharedFirstInSlave ){ mbv[*v].s_rtail = false; mbv[sv].s_ltail = false; }
				else { mbv[*v].s_ltail = false; mbv[sv].s_rtail = false; }

				boost::remove_edge( se, graph );
			}
		}
	}
	
	// discard graphs with bubbles
	if( graph.hasBubbles() ) return false;
	
	// handle unsolvable bifurcations
	for( VertexIterator v = vbegin; v != vend; v++ ) 
	{
		const std::list<Block> &blocks_list = graph.getBlocks(*v);
		const Block &firstBlock = blocks_list.front();
		
		int32_t masterStart = std::min( blocks_list.front().getMasterFrame().getBegin(), blocks_list.back().getMasterFrame().getBegin() );
		int32_t slaveStart = std::min( blocks_list.front().getSlaveFrame().getBegin(), blocks_list.back().getSlaveFrame().getBegin() );
		
		int32_t in_deg = boost::in_degree(*v,graph);
		int32_t out_deg = boost::out_degree(*v,graph);
		
		if( in_deg < 2 && out_deg < 2 ) continue;
		
		Vertex mv, sv, ov; // master/slave linked vertices
		
		if( in_deg == 2 )
		{
			OutEdgeIterator oe_begin, oe_end;
			boost::tie(oe_begin,oe_end) = boost::out_edges(*v,graph);
			if( oe_begin != oe_end ) ov = boost::target(*oe_begin,graph);
			
			InEdgeIterator ebegin,eend, me, se;
			boost::tie(ebegin,eend) = boost::in_edges(*v,graph);
			
			for( InEdgeIterator e = ebegin; e != eend; e++ )
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e);
		
				if( edge_prop.kind == MASTER_EDGE )
				{
					me = e;
					mv = boost::source(*e,graph); 
				}
				else
				{
					se = e;
					sv = boost::source(*e,graph);
				}
			}
			
			const std::list<Block> &masterNextBlocks = graph.getBlocks(mv);
			const std::list<Block> &slaveNextBlocks = graph.getBlocks(sv);
			
			int32_t nextMasterStart = std::min( masterNextBlocks.front().getMasterFrame().getBegin(), masterNextBlocks.back().getMasterFrame().getBegin() );
			int32_t nextSlaveStart = std::min( slaveNextBlocks.front().getSlaveFrame().getBegin(), slaveNextBlocks.back().getSlaveFrame().getBegin() );
			
			bool sharedFirstInMaster = masterStart <= nextMasterStart;
			bool sharedFirstInSlave = slaveStart <= nextSlaveStart;
			
			if( oe_begin != oe_end )
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *oe_begin);
				
				if( edge_prop.kind == MASTER_EDGE )
				{
					mbv[*v].valid = false;
					if( sharedFirstInSlave ) mbv[sv].s_ltail = false; else mbv[sv].s_rtail = false;
					boost::remove_edge( *se, graph );
				}
				else // SLAVE_EDGE
				{
					mbv[*v].valid = false;
					
					if( sharedFirstInSlave ){ mbv[sv].s_ltail = false; mbv[ov].s_rtail = false; } 
					else{ mbv[sv].s_rtail = false; mbv[ov].s_ltail = false; }
					
					boost::remove_edge( mv, *v, graph );
					boost::remove_edge( sv, *v, graph );
				}   
			}
			else
			{
				mbv[*v].valid = false;
				if( sharedFirstInSlave ) mbv[sv].s_ltail = false; else mbv[sv].s_rtail = false;
				boost::remove_edge( *se, graph );
			}
		}
		
		if( out_deg == 2 )
		{
			InEdgeIterator ie_begin, ie_end;
			boost::tie(ie_begin,ie_end) = boost::in_edges(*v,graph);
			if( ie_begin != ie_end ) ov = boost::source(*ie_begin,graph);
			
			OutEdgeIterator ebegin,eend, me, se;
			boost::tie(ebegin,eend) = boost::out_edges(*v,graph);
			
			for( OutEdgeIterator e = ebegin; e != eend; e++ )
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e);
		
				if( edge_prop.kind == MASTER_EDGE )
				{
					me = e;
					mv = boost::target(*e,graph); 
				}
				else
				{
					se = e;
					sv = boost::target(*e,graph);
				}
			}
			
			const std::list<Block> &masterNextBlocks = graph.getBlocks(mv);
			const std::list<Block> &slaveNextBlocks = graph.getBlocks(sv);
			
			int32_t nextMasterStart = std::min( masterNextBlocks.front().getMasterFrame().getBegin(), masterNextBlocks.back().getMasterFrame().getBegin() );
			int32_t nextSlaveStart = std::min( slaveNextBlocks.front().getSlaveFrame().getBegin(), slaveNextBlocks.back().getSlaveFrame().getBegin() );
			
			bool sharedFirstInMaster = masterStart <= nextMasterStart;
			bool sharedFirstInSlave = slaveStart <= nextSlaveStart;
			
			if( ie_begin != ie_end )
			{
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *ie_begin);
				
				if( edge_prop.kind == MASTER_EDGE )
				{
					mbv[*v].valid = false;
					if( sharedFirstInSlave ) mbv[sv].s_ltail = false; else mbv[sv].s_rtail = false;
					boost::remove_edge( se, graph );
				}
				else // SLAVE_EDGE
				{
					mbv[*v].valid = false;
					
					if( sharedFirstInSlave ){ mbv[sv].s_ltail = false; mbv[ov].s_rtail = false; } 
					else{ mbv[sv].s_rtail = false; mbv[ov].s_ltail = false; }
					
					boost::remove_edge( *v, mv, graph );
					boost::remove_edge( *v, sv, graph );
				}   
			}
			else
			{
				mbv[*v].valid = false;
				if( sharedFirstInSlave ) mbv[sv].s_ltail = false; else mbv[sv].s_rtail = false;
				boost::remove_edge( se, graph );
			}
		}
	}
	
	return true;
}


bool PctgBuilder::getMergePaths( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::vector<MergeBlock> &mbv, MergeBlockLists &merge_paths ) const
{
	typedef CompactAssemblyGraph::Vertex Vertex;
	typedef CompactAssemblyGraph::VertexIterator VertexIterator;
	typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;
	
	uint64_t out_degree = boost::out_degree(v,graph);
	uint64_t in_degree = boost::in_degree(v,graph);
	
	if( out_degree >= 2 )
	{
		std::cerr << "[error] Found vertex with output degree >= 2 in fork-solved graph (this should NOT happen!) ==> (" << mbv[v].m_id << "," << mbv[v].s_id << ")" << std::endl;
		return false;
	}
	
	if( in_degree >= 2 )
	{
		std::cerr << "[error] Found vertex with input degree >= 2 in fork-solved graph (this should NOT happen!) ==> (" << mbv[v].m_id << "," << mbv[v].s_id << ")" << std::endl;
		return false;
	}
	
	if( mbv[v].valid ) merge_paths.front().push_back( mbv[v] );
	
	if( out_degree == 0 ) return true;
	
	if( out_degree == 1 )
	{
		OutEdgeIterator ebegin,eend;
		boost::tie(ebegin,eend) = boost::out_edges(v,graph);
		
		Vertex v_cur = boost::source(*ebegin,graph);
		Vertex v_nxt = boost::target(*ebegin,graph);
		
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *ebegin);
		double weight = edge_prop.weight;
		
		//TODO: take care also of master edges, and edges with negative weight (i.e. not enough evidences)
		if( mbv[v_nxt].valid && edge_prop.kind == SLAVE_EDGE && weight >= 0 && weight < 0.3 )
		{
			const std::list<Block> &cur_blocks = graph.getBlocks(v_cur);
			const std::list<Block> &nxt_blocks = graph.getBlocks(v_nxt);
			
			int32_t curSlaveStart = std::min( cur_blocks.front().getSlaveFrame().getBegin(), cur_blocks.back().getSlaveFrame().getBegin() );
			int32_t nxtSlaveStart = std::min( nxt_blocks.front().getSlaveFrame().getBegin(), nxt_blocks.back().getSlaveFrame().getBegin() );
			
			bool curFirstInSlave = curSlaveStart <= nxtSlaveStart;
			
			std::cerr << "[debug] Found (linear) SLAVE mis-assembly in ctg " << mbv[v_cur].s_id << std::endl;
			
			if( curFirstInSlave ){ mbv[v_cur].s_rtail = false; mbv[v_nxt].s_ltail = false; }
			else { mbv[v_cur].s_ltail = false; mbv[v_nxt].s_rtail = false; }
			
			merge_paths.push_front( std::list<MergeBlock>() );
			return this->getMergePaths( graph, v_nxt, mbv, merge_paths );
		}
		
		return this->getMergePaths( graph, v_nxt, mbv, merge_paths );
	}
}


void PctgBuilder::getMergePath( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, MergeBlockLists &merge_paths,
								bool prevMisAssembly, EdgeKindType prevEdgeKind, bool sharedFirst ) const
{
	typedef CompactAssemblyGraph::Vertex Vertex;
	typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;
	typedef CompactAssemblyGraph::InEdgeIterator InEdgeIterator;

	const std::list<Block> &blocks_list = graph.getBlocks(v);
	const Block &firstBlock = blocks_list.front();

	int32_t masterStart = std::min( blocks_list.front().getMasterFrame().getBegin(), blocks_list.back().getMasterFrame().getBegin() );
	int32_t slaveStart = std::min( blocks_list.front().getSlaveFrame().getBegin(), blocks_list.back().getSlaveFrame().getBegin() );

	uint64_t out_degree = boost::out_degree(v,graph);
	uint64_t in_degree = boost::in_degree(v,graph);

	MergeBlock merge_block;

	merge_block.vertex = v;
	merge_block.m_id = firstBlock.getMasterId();
	merge_block.s_id = firstBlock.getSlaveId();

	merge_block.ext_slave_next = true;
	merge_block.ext_slave_prev = true;

	merge_block.m_ltail = true;
	merge_block.m_rtail = true;
	merge_block.s_ltail = true;
	merge_block.s_rtail = true;

	if( prevMisAssembly )
	{
		if( prevEdgeKind == MASTER_EDGE )
		{
			if( sharedFirst ) merge_block.m_ltail = false; else merge_block.m_rtail = false;
		}
		else // prevEdgeKind == SLAVE_EDGE
		{
			if( sharedFirst ) merge_block.s_ltail = false; else merge_block.s_rtail = false;
		}
	}


	if( out_degree == 0 ) return merge_paths.front().push_back( merge_block );

	if( out_degree == 1 )
	{
		OutEdgeIterator ebegin,eend;
		boost::tie(ebegin,eend) = boost::out_edges(v,graph);

		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *ebegin);
		Vertex v_next = boost::target(*ebegin,graph);

		merge_paths.front().push_back( merge_block );
		return this->getMergePath( graph, v_next, merge_paths );
	}

	if( out_degree == 2 )
	{
		ForkType fork_type = this->getOutForkType(graph,v);

		Vertex mv, sv;
		OutEdgeIterator ebegin,eend;
		boost::tie(ebegin,eend) = boost::out_edges(v,graph);

		for( OutEdgeIterator e = ebegin; e != eend; e++ )
		{
			EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e);
			if( edge_prop.kind == MASTER_EDGE ) mv = boost::target(*e,graph); else sv = boost::target(*e,graph);
		}

		const std::list<Block> &masterNextBlocks = graph.getBlocks(mv);
		const std::list<Block> &slaveNextBlocks = graph.getBlocks(sv);

		int32_t nextMasterStart = std::min( masterNextBlocks.front().getMasterFrame().getBegin(), masterNextBlocks.back().getMasterFrame().getBegin() );
		int32_t nextSlaveStart = std::min( slaveNextBlocks.front().getSlaveFrame().getBegin(), slaveNextBlocks.back().getSlaveFrame().getBegin() );

		bool sharedFirstInMaster = masterStart <= nextMasterStart;
		bool sharedFirstInSlave = slaveStart <= nextSlaveStart;

		if( fork_type == MIS_MASTER )
		{
			std::cerr << "found MASTER misassembly in ctg "  << firstBlock.getMasterId() << std::endl;

			if( sharedFirstInMaster ) merge_block.m_rtail = false; else merge_block.m_ltail = false;
			//merge_block.out_m_mis = true;

			merge_paths.front().push_back( merge_block );
			this->getMergePath( graph, sv, merge_paths );

			merge_paths.push_front( std::list<MergeBlock>() );
			return this->getMergePath( graph, mv, merge_paths, true, MASTER_EDGE, sharedFirstInMaster );
		}
		else if( fork_type == MIS_SLAVE )
		{
			std::cerr << "found SLAVE misassembly in ctg " << firstBlock.getSlaveId() << std::endl;

			if( sharedFirstInSlave ) merge_block.s_rtail = false; else merge_block.s_ltail = false;
			//merge_block.out_s_mis = true;

			merge_paths.front().push_back( merge_block );
			this->getMergePath( graph, mv, merge_paths );

			merge_paths.push_front( std::list<MergeBlock>() );
			return this->getMergePath( graph, sv, merge_paths, true, SLAVE_EDGE, sharedFirstInSlave );
		}
		else // fork_type == UNKNOWN
		{
			std::list< MergeBlock > &path = merge_paths.front();

			InEdgeIterator ie_it,ie_end;

			boost::tie(ie_it,ie_end) = boost::in_edges(v,graph);
			while( ie_it != ie_end )
			{
				const std::list<Block> &bl = graph.getBlocks( boost::source(*ie_it,graph) );
				const Block &b = blocks_list.front();

				if( merge_block.m_id != b.getMasterId() && merge_block.s_id != b.getSlaveId() )
				{
					path.back().ext_slave_next = false;
					break;
				}

				path.pop_back();
				boost::tie(ie_it,ie_end) = boost::in_edges( boost::source(*ie_it,graph), graph );
			}

			OutEdgeIterator oe_it,oe_end;

			boost::tie(oe_it,oe_end) = boost::out_edges(sv,graph);
			while( oe_it != oe_end )
			{
				sv = boost::target(*oe_it,graph);
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *oe_it);

				if( edge_prop.kind == MASTER_EDGE )
				{
					merge_paths.push_front( std::list<MergeBlock>() );
					this->getMergePath( graph, sv, merge_paths );
					break;
				}

				boost::tie(oe_it,oe_end) = boost::out_edges(sv,graph);
			}

			merge_paths.push_front( std::list<MergeBlock>() );
			return this->getMergePath( graph, mv, merge_paths );
		}
	} // end if out_degree == 2
}

/* void PctgBuilder::getMergeLists( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex v, MergeDescriptorLists &mergeLists,
								 bool hasMisAssembly, EdgeKindType edgeKind, bool sharedFirst ) const
{
    typedef CompactAssemblyGraph::Vertex Vertex;
	typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;

    const std::list<Block> &blocks_list = graph.getBlocks(v);

    const Block &firstBlock = blocks_list.front();
	const Block &lastBlock = blocks_list.back();

	const Frame &firstMasterFrame = firstBlock.getMasterFrame();
	const Frame &firstSlaveFrame = firstBlock.getSlaveFrame();
	const Frame &lastMasterFrame = lastBlock.getMasterFrame();
	const Frame &lastSlaveFrame = lastBlock.getSlaveFrame();

    int32_t masterId = firstMasterFrame.getContigId();
    int32_t slaveId = firstSlaveFrame.getContigId();

    int32_t masterStart = std::min( firstMasterFrame.getBegin(), lastMasterFrame.getBegin() );
    int32_t masterEnd = std::max( firstMasterFrame.getEnd(), lastMasterFrame.getEnd() );
    int32_t slaveStart = std::min( firstSlaveFrame.getBegin(), lastSlaveFrame.getBegin() );
    int32_t slaveEnd = std::max( firstSlaveFrame.getEnd(), lastSlaveFrame.getEnd() );

    // load contig that should be merged with pctg.
    Contig masterCtg = this->loadMasterContig(masterId);
    Contig slaveCtg = this->loadSlaveContig(slaveId);

	uint64_t out_degree = boost::out_degree(v,graph);

	// find best alignment between the contigs
	BestCtgAlignment *bestAlign = new BestCtgAlignment();
	this->findBestAlignment( *bestAlign, masterCtg, masterStart, masterEnd, slaveCtg, slaveStart, slaveEnd, blocks_list );

	MergeDescriptor mergeDesc;
	mergeDesc.align = bestAlign;

	// find start/end positions of the best alignment
	std::pair<uint64_t,uint64_t> alignStart, alignEnd, alignStartTmp, alignEndTmp;

	if( bestAlign->main_homology() >= MIN_HOMOLOGY )
	{
		first_match_pos( bestAlign->at(0), alignStart );
		last_match_pos( bestAlign->at(bestAlign->size()-1), alignEnd );

		// compute masterCtg tails (i1,i2) and slaveCtg tails (j1,j2)
		uint64_t i1 = alignStart.first;
		uint64_t i2 = masterCtg.size() - alignEnd.first - 1;
		uint64_t j1 = alignStart.second;
		uint64_t j2 = slaveCtg.size() - alignEnd.second - 1;

		const MyAlignment& left = bestAlign->left();
		const MyAlignment& right = bestAlign->right();

		size_t mt = 0.3 * masterCtg.size();
		size_t st = 0.3 * slaveCtg.size();

		uint64_t left_min_len = 0.7 * std::min(i1,j1);
		uint64_t right_min_len = 0.7 * std::min(i2,j2);
		uint64_t threshold = std::min( size_t(500), std::min(mt,st) );

		if( std::min(i1,j1) >= threshold && this->is_good(left,left_min_len) )
		{
			first_match_pos(left,alignStart);
			if( bestAlign->is_left_rev() ) std::swap( alignStart.first, alignStart.second );
		}

		if( std::min(i2,j2) >= threshold && this->is_good(right,right_min_len) )
		{
			alignEndTmp = alignEnd;
			last_match_pos(right,alignEnd);

			if( bestAlign->is_right_rev() ) std::swap( alignEnd.first, alignEnd.second );

			alignEnd.first += alignEndTmp.first+1;
			alignEnd.second += alignEndTmp.second+1;
		}
	}
	else // bad alignment between blocks
	{
		alignStart = std::make_pair(masterStart,slaveStart);
		alignEnd = std::make_pair(masterEnd,slaveEnd);
	}

	if( bestAlign->isCtgReversed() )
	{
		uint64_t tempPos = alignStart.second;
		alignStart.second = slaveCtg.size() - alignEnd.second - 1;
		alignEnd.second = slaveCtg.size() - tempPos - 1;
	}

	mergeDesc.masterMainAlign.first = alignStart.first;
	mergeDesc.masterMainAlign.second = alignEnd.first;
	mergeDesc.slaveMainAlign.first = alignStart.second;
	mergeDesc.slaveMainAlign.second = alignEnd.second;

	if( out_degree == 2 ) // handle putative misassembly
	{
		ForkType type = this->getForkType(graph,v);

		mergeDesc.mergeType = FORK_MERGE;
		mergeDesc.forkType = type;

		Vertex mv, sv; // next master/slave linked vertices

		OutEdgeIterator ebegin,eend;
		boost::tie(ebegin,eend) = boost::out_edges(v,graph);

		for( OutEdgeIterator e = ebegin; e != eend; e++ )
		{
			EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e);
			if( edge_prop.kind == MASTER_EDGE ) mv = boost::target(*e,graph); else sv = boost::target(*e,graph);
		}

		if( type == MIS_MASTER || type == MIS_SLAVE )
		{
			const std::list<Block> &masterNextBlocks = graph.getBlocks(mv);
			const std::list<Block> &slaveNextBlocks = graph.getBlocks(sv);

			int32_t nextMasterStart = std::min( masterNextBlocks.front().getMasterFrame().getBegin(), masterNextBlocks.back().getMasterFrame().getBegin() );
			int32_t nextSlaveStart = std::min( slaveNextBlocks.front().getSlaveFrame().getBegin(), slaveNextBlocks.back().getSlaveFrame().getBegin() );

			bool sharedFirstInMaster = masterStart <= nextMasterStart;
			bool sharedFirstInSlave = slaveStart <= nextSlaveStart;

			if( type == MIS_MASTER )
			{
				std::cerr << "found MASTER misassembly "  << firstMasterFrame.getContigId() << std::endl;

				mergeDesc.masterPos.first  = sharedFirstInMaster ? 0 : alignStart.first;
				mergeDesc.masterPos.second = sharedFirstInMaster ? alignEnd.first : masterCtg.size()-1;
				mergeLists.front().push_back(mergeDesc);

				this->getMergeLists( graph, sv, mergeLists, false );

				mergeLists.push_front( std::list<MergeDescriptor>() );
				this->getMergeLists( graph, mv, mergeLists, true, MASTER_EDGE, sharedFirstInMaster );
			}
			else // type == MIS_SLAVE
			{
				std::cerr << "found SLAVE misassembly" << std::endl;

				mergeDesc.slavePos.first  = sharedFirstInSlave ? 0 : alignStart.second;
				mergeDesc.slavePos.second = sharedFirstInSlave ? alignEnd.second : slaveCtg.size()-1;
				mergeLists.front().push_back(mergeDesc);

				this->getMergeLists( graph, mv, mergeLists, false );

				mergeLists.push_front( std::list<MergeDescriptor>() );
				this->getMergeLists( graph, sv, mergeLists, true, SLAVE_EDGE, sharedFirstInSlave );
			}
		}
		else // type == UNKNOWN || type == REPEAT
		{
			mergeLists.front().push_back(mergeDesc);

			OutEdgeIterator e_it,e_end;

			boost::tie(e_it,e_end) = boost::out_edges(mv,graph);
			while( e_it != e_end )
			{
				mv = boost::target(*e_it,graph);
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e_it); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e_it);

				if( edge_prop.kind == MASTER_EDGE )
				{
					mergeLists.push_front( std::list<MergeDescriptor>() );
					this->getMergeLists( graph, mv, mergeLists, false );
					break;
				}

				boost::tie(e_it,e_end) = boost::out_edges(mv,graph);
			}

			boost::tie(e_it,e_end) = boost::out_edges(sv,graph);
			while( e_it != e_end )
			{
				sv = boost::target(*e_it,graph);
				EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e_it); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e_it);

				if( edge_prop.kind == MASTER_EDGE )
				{
					mergeLists.push_front( std::list<MergeDescriptor>() );
					this->getMergeLists( graph, sv, mergeLists, false );
					break;
				}

				boost::tie(e_it,e_end) = boost::out_edges(sv,graph);
			}


			return;
		}
	}
	else // out_degree == 0 || out_degree == 1
	{
		mergeDesc.mergeType = LINEAR_MERGE;

		if( hasMisAssembly )
		{
			if( edgeKind == MASTER_EDGE )
			{
				mergeDesc.masterPos.first  = sharedFirst ? alignStart.first : 0;
				mergeDesc.masterPos.second = sharedFirst ? masterCtg.size()-1 : alignEnd.first;
			}
			else
			{
				mergeDesc.slavePos.first  = sharedFirst ? alignStart.second : 0;
				mergeDesc.slavePos.second = sharedFirst ? slaveCtg.size()-1 : alignEnd.second;
			}
		}

		mergeLists.front().push_back(mergeDesc);

		OutEdgeIterator ebegin,eend;
		boost::tie(ebegin,eend) = boost::out_edges(v,graph);

		for( OutEdgeIterator e = ebegin; e != eend; e++ )
		{
			EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e);
			this->getMergeLists( graph, boost::target(*e,graph), mergeLists, false );
		}
	}
}

*/

ForkType PctgBuilder::getInForkType( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v ) const
{
	typedef CompactAssemblyGraph::InEdgeIterator InEdgeIterator;
	typedef CompactAssemblyGraph::Vertex Vertex;

	Vertex mv, sv; // master/slave linked vertices

	InEdgeIterator ebegin,eend;
	boost::tie(ebegin,eend) = boost::in_edges(v,graph);

	for( InEdgeIterator e = ebegin; e != eend; e++ )
	{
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e);
		if( edge_prop.kind == MASTER_EDGE ) mv = boost::source(*e,graph); else sv = boost::source(*e,graph);
	}

	return getForkType( graph, v, mv, sv );
}


ForkType PctgBuilder::getOutForkType( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v ) const
{
	typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;
	typedef CompactAssemblyGraph::Vertex Vertex;

	Vertex mv, sv; // master/slave linked vertices

	OutEdgeIterator ebegin,eend;
	boost::tie(ebegin,eend) = boost::out_edges(v,graph);

	for( OutEdgeIterator e = ebegin; e != eend; e++ )
	{
		EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *e); //EdgeKindType kind = boost::get(boost::edge_kind_t(), graph, *e);
		if( edge_prop.kind == MASTER_EDGE ) mv = boost::target(*e,graph); else sv = boost::target(*e,graph);
	}

	return getForkType( graph, v, mv, sv );
}


ForkType PctgBuilder::getForkType( const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v,
								   CompactAssemblyGraph::Vertex &mv, CompactAssemblyGraph::Vertex &sv ) const
{
	const std::list<Block>& sharedBlocks = graph.getBlocks(v);
	const std::list<Block>& masterBlocks = graph.getBlocks(mv);
	const std::list<Block>& slaveBlocks = graph.getBlocks(sv);

	int32_t nh, xt;

	std::vector< uint64_t > m_pairs( masterBam.size(), 0 );
	std::vector< uint64_t > m_cmp( masterBam.size(), 0 );
	std::vector< uint64_t > m_pairs_unmap( masterBam.size(), 0 );
	std::vector< uint64_t > m_pairs_ormap( masterBam.size(), 0 );

	std::vector< uint64_t > m_mates( masterMpBam.size(), 0 );
	std::vector< uint64_t > m_cmm( masterMpBam.size(), 0 );
	std::vector< uint64_t > m_mates_unmap( masterMpBam.size(), 0 );
	std::vector< uint64_t > m_mates_ormap( masterMpBam.size(), 0 );
    
    //std::vector< uint64_t > m_mates_cmap( masterMpBam.size(), 0 );
    //std::vector< uint64_t > m_mates_wmap( masterMpBam.size(), 0 );

	std::vector< uint64_t > s_pairs( slaveBam.size(), 0 );
	std::vector< uint64_t > s_cmp( slaveBam.size(), 0 );
	std::vector< uint64_t > s_pairs_unmap( slaveBam.size(), 0 );
	std::vector< uint64_t > s_pairs_ormap( slaveBam.size(), 0 );

	std::vector< uint64_t > s_mates( slaveMpBam.size(), 0 );
	std::vector< uint64_t > s_cmm( slaveMpBam.size(), 0 );
	std::vector< uint64_t > s_mates_unmap( slaveMpBam.size(), 0 );
	std::vector< uint64_t > s_mates_ormap( slaveMpBam.size(), 0 );
    
    //std::vector< uint64_t > s_mates_cmap( slaveMpBam.size(), 0 );
    //std::vector< uint64_t > s_mates_wmap( slaveMpBam.size(), 0 );

	const int32_t masterId = sharedBlocks.front().getMasterId();
	const int32_t slaveId = sharedBlocks.front().getSlaveId();

	const int32_t sharedMasterStart = std::min( sharedBlocks.front().getMasterFrame().getBegin(), sharedBlocks.back().getMasterFrame().getBegin() );
	const int32_t sharedMasterEnd = std::max( sharedBlocks.front().getMasterFrame().getEnd(), sharedBlocks.back().getMasterFrame().getEnd() );
	const int32_t sharedSlaveStart = std::min( sharedBlocks.front().getSlaveFrame().getBegin(), sharedBlocks.back().getSlaveFrame().getBegin() );
	const int32_t sharedSlaveEnd = std::max( sharedBlocks.front().getSlaveFrame().getEnd(), sharedBlocks.back().getSlaveFrame().getEnd() );

	const int32_t linkedMasterStart = std::min( masterBlocks.front().getMasterFrame().getBegin(), masterBlocks.back().getMasterFrame().getBegin() );
	const int32_t linkedMasterEnd = std::min( masterBlocks.front().getMasterFrame().getEnd(), masterBlocks.back().getMasterFrame().getEnd() );
	const int32_t linkedSlaveStart = std::min( slaveBlocks.front().getSlaveFrame().getBegin(), slaveBlocks.back().getSlaveFrame().getBegin() );
	const int32_t linkedSlaveEnd = std::min( slaveBlocks.front().getSlaveFrame().getEnd(), slaveBlocks.back().getSlaveFrame().getEnd() );

	bool sharedFirstInMaster = sharedMasterStart <= linkedMasterStart;
	bool sharedFirstInSlave = sharedSlaveStart <= linkedSlaveStart;

	// limit on master/slave from where to search for the pairs/mates
	int32_t k_m = sharedFirstInMaster ? sharedMasterEnd+1 : sharedMasterStart-1;
	int32_t k_s = sharedFirstInSlave ? sharedSlaveEnd+1 : sharedSlaveStart-1;

	uint64_t r1_m = sharedFirstInMaster ? sharedMasterEnd + 1 : _masterRef->at(masterId).RefLength - sharedMasterStart;
	uint64_t r1_s = sharedFirstInSlave ? sharedSlaveEnd + 1 : _slaveRef->at(slaveId).RefLength - sharedSlaveStart;
	uint64_t r2_m = sharedFirstInMaster ? _masterRef->at(masterId).RefLength - sharedMasterEnd - 1 : sharedMasterStart;
	uint64_t r2_s = sharedFirstInSlave ? _slaveRef->at(slaveId).RefLength - sharedSlaveEnd - 1 : sharedSlaveStart;

	uint64_t r1_min = r1_m <= r1_s ? r1_m : r1_s;
	uint64_t r2_min = r2_m <= r2_s ? r2_m : r2_s;

	// gap between the blocks in master and slave
	int32_t masterGap = sharedFirstInMaster ? (linkedMasterStart - sharedMasterEnd) : (sharedMasterStart - linkedMasterEnd);
	int32_t slaveGap = sharedFirstInSlave ? (linkedSlaveStart - sharedSlaveEnd) : (sharedSlaveStart - linkedSlaveEnd);

	// COMPUTE STATISTICS FOR PE LIBRARIES
	for( int i=0; i < masterBam.size(); i++ )
	{
		BamAlignment align;

		// retrieve BAM readers for master/slave for the current PE library
		BamReader *masterBamReader = masterBam.getBamReader(i);
		BamReader *slaveBamReader = slaveBam.getBamReader(i);

		double masterLibMean = masterBam.getISizeMean(i);
		double masterLibStd = masterBam.getISizeStd(i);
		double slaveLibMean = slaveBam.getISizeMean(i);
		double slaveLibStd = slaveBam.getISizeStd(i);

		uint64_t r1 = r1_min;

		if( masterLibMean - 2*masterLibStd > r1+r2_min ) continue;
		if( slaveLibMean - 2*slaveLibStd > r1+r2_min ) continue;

		//if( masterLibMean - 2*masterLibStd > r1+r2_min ) std::cerr << "master PE stats may be unreliable for lib " << masterLibMean << "(±" << masterLibStd << ") r1+r2=" << r1+r2_min << "bp" << std::endl;
		//if( slaveLibMean - 2*slaveLibStd > r1+r2_min ) std::cerr << "slave PE stats may be unreliable for lib " << slaveLibMean << "(±" << slaveLibStd << ") r1+r2=" << r1+r2_min << "bp" << std::endl;

		// reduce interval if longer than insert size (mean + 2*std)
		if( masterLibMean + 2*masterLibStd < r1 ) r1 = masterLibMean + 2*masterLibStd;
		if( slaveLibMean + 2*slaveLibStd < r1 ) r1 = slaveLibMean + 2*masterLibStd;

		int32_t masterStart = sharedFirstInMaster ? (k_m - r1) : (k_m + 1);
		int32_t masterEnd = sharedFirstInMaster ? k_m : (k_m + r1 + 1);
		int32_t slaveStart = sharedFirstInSlave ? (k_s - r1) : (k_s + 1);
		int32_t slaveEnd = sharedFirstInSlave ? k_s : (k_s + r1 + 1);

		// lock current master BAM file
		masterBam.lockBamReader(i);

		// set region in master BAM reader
		masterBamReader->SetRegion( masterId, masterStart, masterId, masterEnd );
		// compute pair-end statistics for the master
		while( masterBamReader->GetNextAlignmentCore(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

			align.BuildCharData(); // fill string fields

			// se la molteplicità non è stata definita, assumo che sia pari ad 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
			if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;

			int32_t maxInsert = masterLibMean + 2*masterLibStd;
			int32_t minInsert = masterLibMean - 2*masterLibStd;

			if( !align.IsMateMapped() )
			{
				// I want to count only reads whose mates are supposed to be mapped after
				if( sharedFirstInMaster && !align.IsReverseStrand() ) m_pairs_unmap[i]++;
				if( !sharedFirstInMaster && align.IsReverseStrand() ) m_pairs_unmap[i]++;
			}

			// discard reads whose mate should be reasonably mapped in another contig
			if( align.IsMateMapped() && align.RefID != align.MateRefID )
			{
				if( sharedFirstInMaster && !align.IsReverseStrand() )
				{
					int32_t endPos = align.Position + minInsert;
					int32_t endPos2 = align.Position + maxInsert;
					int32_t ctgLen = _masterRef->at(masterId).RefLength;
					if( endPos >= ctgLen || endPos2 >= ctgLen ) continue; else m_pairs_ormap[i]++;
				}

				if( !sharedFirstInMaster && align.IsReverseStrand() )
				{
					int32_t startPos = align.Position + readLength - minInsert;
					int32_t startPos2 = align.Position + readLength - maxInsert;
					if( startPos < 0 || startPos2 < 0 ) continue; else m_pairs_ormap[i]++;
				}
			}

			// discard unpaired reads or reads not inside the region
			if( !align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

			if( sharedFirstInMaster )
			{
				if( startRead < masterStart ) continue;

				// read completamente incluse nell'intervallo
				if( endRead < masterEnd )
				{
					if( align.MatePosition >= masterEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() ) m_pairs[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition >= masterEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) m_cmp[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate è mappata completamente oltre l'intervallo
				if( align.MatePosition >= masterEnd )
				{
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() ) m_pairs[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) m_cmp[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate si sovrappone anch'essa all'intervallo
				if( startRead > align.MatePosition )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) m_pairs[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) m_cmp[i]++;

					continue;
				}
			}
			else
			{
				if( endRead >= masterEnd ) continue;

				// read completamente incluse nell'intervallo
				if( startRead >= masterStart )
				{
					if( align.MatePosition < masterStart && align.IsReverseStrand() && !align.IsMateReverseStrand() ) m_pairs[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition < masterStart && align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) m_cmp[i]++;

					continue;
				}

				// read che escono dall'intervallo a sinistra e la cui mate si trova a sinistra
				if( align.MatePosition < startRead )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) m_pairs[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) m_cmp[i]++;

					continue;
				}
			}
		}

		// unlock current master BAM file
		masterBam.unlockBamReader(i);

		// lock current slave BAM file
		slaveBam.lockBamReader(i);

		// set region in slave BAM reader
		slaveBamReader->SetRegion( slaveId, slaveStart, slaveId, slaveEnd );
		// compute pair-end statistics for the slave
		while( slaveBamReader->GetNextAlignmentCore(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

			align.BuildCharData(); // fill string fields

			// se la molteplicità non è stata definita, assumo che sia pari ad 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
			if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;

			int32_t maxInsert = slaveLibMean + 2*slaveLibStd;
			int32_t minInsert = slaveLibMean - 2*slaveLibStd;

			if( !align.IsMateMapped() )
			{
				// I want to count only reads whose mates are supposed to be mapped after
				if( sharedFirstInSlave && !align.IsReverseStrand() ) s_pairs_unmap[i]++;
				if( !sharedFirstInSlave && align.IsReverseStrand() ) s_pairs_unmap[i]++;
			}

			// discard reads whose mate should be mapped in another contig
			if( align.IsMateMapped() && align.RefID != align.MateRefID )
			{
				if( sharedFirstInSlave && !align.IsReverseStrand() )
				{
					int32_t endPos = align.Position + minInsert;
					int32_t endPos2 = align.Position + maxInsert;
					int32_t ctgLen = _slaveRef->at(slaveId).RefLength;
					if( endPos >= ctgLen || endPos2 >= ctgLen ) continue; else s_pairs_ormap[i]++;
				}

				if( !sharedFirstInSlave && align.IsReverseStrand() )
				{
					int32_t startPos = align.Position + readLength - minInsert;
					int32_t startPos2 = align.Position + readLength - maxInsert;
					if( startPos < 0 || startPos2 < 0 ) continue; else s_pairs_ormap[i]++;
				}
			}

			// discard unpaired reads or reads not inside the region
			if( !align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

			if( sharedFirstInSlave )
			{
				if( startRead < slaveStart ) continue;

				// read completamente incluse nell'intervallo
				if( endRead < slaveEnd )
				{
					if( align.MatePosition >= slaveEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() ) s_pairs[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition >= slaveEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) s_cmp[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate è mappata completamente oltre l'intervallo
				if( align.MatePosition >= slaveEnd )
				{
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() ) s_pairs[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) s_cmp[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate si sovrappone anch'essa all'intervallo
				if( startRead > align.MatePosition )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) s_pairs[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) s_cmp[i]++;

					continue;
				}
			}
			else
			{
				if( endRead >= slaveEnd ) continue;

				// read completamente incluse nell'intervallo
				if( startRead >= slaveStart )
				{
					if( align.MatePosition < slaveStart && align.IsReverseStrand() && !align.IsMateReverseStrand() ) s_pairs[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition < slaveStart && align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) s_cmp[i]++;

					continue;
				}

				// read che escono dall'intervallo a sinistra e la cui mate si trova a sinistra
				if( align.MatePosition < startRead )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) s_pairs[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) s_cmp[i]++;

					continue;
				}
			}
		}

		// unlock current slave BAM file
		slaveBam.unlockBamReader(i);

	}

	// COMPUTE STATISTICS FOR MP LIBRARIES
	for( int i=0; i < masterMpBam.size(); i++ )
	{
		BamAlignment align;

		// retrieve BAM readers for master/slave for the current MP library
		BamReader *masterBamReader = masterMpBam.getBamReader(i);
		BamReader *slaveBamReader = slaveMpBam.getBamReader(i);

		double masterLibMean = masterMpBam.getISizeMean(i);
		double masterLibStd = masterMpBam.getISizeStd(i);
		double slaveLibMean = slaveMpBam.getISizeMean(i);
		double slaveLibStd = slaveMpBam.getISizeStd(i);

		uint64_t r1 = r1_min;

		if( masterLibMean - 2*masterLibStd > r1+r2_min ) continue;
		if( slaveLibMean - 2*slaveLibStd > r1+r2_min ) continue;
		//if( masterLibMean - 2*masterLibStd > r1+r2_min ) std::cerr << "master MP stats may be unreliable for lib " << masterLibMean << "(±" << masterLibStd << ") r1+r2=" << r1+r2_min << "bp" << std::endl;
		//if( slaveLibMean - 2*slaveLibStd > r1+r2_min ) std::cerr << "slave MP stats may be unreliable for lib " << slaveLibMean << "(±" << slaveLibStd << ") r1+r2=" << r1+r2_min << "bp" << std::endl;

		if( masterLibMean + 2*masterLibStd < r1 ) r1 = masterLibMean + 2*masterLibStd;
		if( slaveLibMean + 2*slaveLibStd < r1 ) r1 = slaveLibMean + 2*slaveLibStd;

		int32_t masterStart = sharedFirstInMaster ? (k_m - r1) : (k_m + 1);
		int32_t masterEnd = sharedFirstInMaster ? k_m : (k_m + r1 + 1);
		int32_t slaveStart = sharedFirstInSlave ? (k_s - r1) : (k_s + 1);
		int32_t slaveEnd = sharedFirstInSlave ? k_s : (k_s + r1 + 1);

		// lock current master BAM file
		masterMpBam.lockBamReader(i);

		// set region in master BAM reader
		masterBamReader->SetRegion( masterId, masterStart, masterId, masterEnd );
		// compute mate-pairs statistics for the master
		while( masterBamReader->GetNextAlignmentCore(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

			align.BuildCharData(); // fill string fields

			// se la molteplicità non è stata definita, assumo che sia pari ad 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
			if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;

			int32_t maxInsert = masterLibMean + 2*masterLibStd;
			int32_t minInsert = masterLibMean - 2*masterLibStd;

			if( !align.IsMateMapped() )
			{
				// I want to count only reads whose mates are supposed to be mapped after
				if( sharedFirstInMaster && !align.IsReverseStrand() ) m_mates_unmap[i]++;
				if( !sharedFirstInMaster && align.IsReverseStrand() ) m_mates_unmap[i]++;
			}

			// discard reads whose mate should be reasonably mapped in another contig
			if( align.IsMateMapped() && align.RefID != align.MateRefID )
			{
				if( sharedFirstInMaster && !align.IsReverseStrand() )
				{
					int32_t endPos = align.Position + minInsert;
					int32_t endPos2 = align.Position + maxInsert;
					int32_t ctgLen = _masterRef->at(masterId).RefLength;
					if( endPos >= ctgLen || endPos2 >= ctgLen ) continue; else m_mates_ormap[i]++;
				}

				if( !sharedFirstInMaster && align.IsReverseStrand() )
				{
					int32_t startPos = align.Position + readLength - minInsert;
					int32_t startPos2 = align.Position + readLength - maxInsert;
					if( startPos < 0 || startPos2 < 0 ) continue; else m_mates_ormap[i]++;
				}
			}

			// discard unpaired reads or reads not inside the region
			if( !align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

			if( sharedFirstInMaster )
			{
				if( startRead < masterStart ) continue;

				// read completamente incluse nell'intervallo
				if( endRead < masterEnd )
				{
					if( align.MatePosition >= masterEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() ) m_mates[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition >= masterEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) m_cmm[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate è mappata completamente oltre l'intervallo
				if( align.MatePosition >= masterEnd )
				{
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() ) m_mates[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) m_cmm[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate si sovrappone anch'essa all'intervallo
				if( startRead > align.MatePosition )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) m_mates[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) m_cmm[i]++;

					continue;
				}
			}
			else
			{
				if( endRead >= masterEnd ) continue;

				// read completamente incluse nell'intervallo
				if( startRead >= masterStart )
				{
					if( align.MatePosition < masterStart && align.IsReverseStrand() && !align.IsMateReverseStrand() ) m_mates[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition < masterStart && align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) m_cmm[i]++;

					continue;
				}

				// read che escono dall'intervallo a sinistra e la cui mate si trova a sinistra
				if( align.MatePosition < startRead )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) m_mates[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) m_cmm[i]++;

					continue;
				}
			}
		}

		masterMpBam.unlockBamReader(i);
		slaveMpBam.lockBamReader(i);

		// set region in slave BAM reader
		slaveBamReader->SetRegion( slaveId, slaveStart, slaveId, slaveEnd );
		// compute mate-pairs statistics for the slave
		while( slaveBamReader->GetNextAlignmentCore(align) )
		{
			// discard bad quality reads
			if( !align.IsMapped() || !align.IsPaired() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

			align.BuildCharData(); // fill string fields

			// se la molteplicità non è stata definita, assumo che sia pari ad 1
			if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
			if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
			if( nh != 1 || xt != 'U' ) continue; // discard reads with multiplicity greater than 1

			int32_t readLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t endRead = startRead + readLength - 1;

			int32_t maxInsert = slaveLibMean + 2*slaveLibStd;
			int32_t minInsert = slaveLibMean - 2*slaveLibStd;

			if( !align.IsMateMapped() )
			{
				// I want to count only reads whose mates are supposed to be mapped after
				if( sharedFirstInSlave && !align.IsReverseStrand() ) s_mates_unmap[i]++;
				if( !sharedFirstInSlave && align.IsReverseStrand() ) s_mates_unmap[i]++;
			}

			// discard reads whose mate should be mapped in another contig
			if( align.IsMateMapped() && align.RefID != align.MateRefID )
			{
				if( sharedFirstInSlave && !align.IsReverseStrand() )
				{
					int32_t endPos = align.Position + minInsert;
					int32_t endPos2 = align.Position + maxInsert;
					int32_t ctgLen = _slaveRef->at(slaveId).RefLength;
					if( endPos >= ctgLen || endPos2 >= ctgLen ) continue; else s_mates_ormap[i]++;
				}

				if( !sharedFirstInSlave && align.IsReverseStrand() )
				{
					int32_t startPos = align.Position + readLength - minInsert;
					int32_t startPos2 = align.Position + readLength - maxInsert;
					if( startPos < 0 || startPos2 < 0 ) continue; else s_mates_ormap[i]++;
				}
			}

			// discard unpaired reads or reads not inside the region
			if( !align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

			if( sharedFirstInSlave )
			{
				if( startRead < slaveStart ) continue;

				// read completamente incluse nell'intervallo
				if( endRead < slaveEnd )
				{
					if( align.MatePosition >= slaveEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() ) s_mates[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition >= slaveEnd && !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) s_cmm[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate è mappata completamente oltre l'intervallo
				if( align.MatePosition >= slaveEnd )
				{
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() ) s_mates[i]++;

					int32_t isize = align.MatePosition + readLength - startRead;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() && good_isize ) s_cmm[i]++;

					continue;
				}

				// read che escono dall'intervallo a destra e la cui mate si sovrappone anch'essa all'intervallo
				if( startRead > align.MatePosition )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) s_mates[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) s_cmm[i]++;

					continue;
				}
			}
			else
			{
				if( endRead >= slaveEnd ) continue;

				// read completamente incluse nell'intervallo
				if( startRead >= slaveStart )
				{
					if( align.MatePosition < slaveStart && align.IsReverseStrand() && !align.IsMateReverseStrand() ) s_mates[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.MatePosition < slaveStart && align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) s_cmm[i]++;

					continue;
				}

				// read che escono dall'intervallo a sinistra e la cui mate si trova a sinistra
				if( align.MatePosition < startRead )
				{
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() ) s_mates[i]++;

					int32_t isize = endRead - align.MatePosition + 1;
					bool good_isize = (isize <= maxInsert) && (isize >= minInsert);
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() && good_isize ) s_cmm[i]++;

					continue;
				}
			}
		}

		slaveMpBam.unlockBamReader(i);
	}

	uint32_t masterMisEvid = 0, slaveMisEvid = 0;

	for( int i=0; i < m_pairs.size(); i++ )
	{
		uint64_t masterTotPairs = m_pairs[i] + m_pairs_unmap[i] + m_pairs_ormap[i];
		uint64_t slaveTotPairs = s_pairs[i] + s_pairs_unmap[i] + s_pairs_ormap[i];
		
		if(slaveTotPairs >= 10 || masterTotPairs >= 10)
		{
			double masterScore = (masterTotPairs > 0) ? 100 * m_pairs[i] / double(masterTotPairs) : 0.0;
			double slaveScore = (slaveTotPairs > 0) ? 100 * s_pairs[i] / double(slaveTotPairs) : 0.0;
			
			if( masterScore >= 95.0 && slaveScore <= 5.0 ) slaveMisEvid++;
			if( slaveScore >= 95.0 && masterScore <= 5.0 ) masterMisEvid++;
		}
	}
    
    /*std::cerr << "MISASSEMBLY IN (" << masterId << "," << slaveId << ") PAIRS-STATS" << std::endl;
    std::cerr << "Master: paired=" << m_pairs[0] << " unpaired=" << m_pairs_unmap[0] << " paired_other_ctg=" << m_pairs_ormap[0] << std::endl;
    std::cerr << "Slave:  paired=" << s_pairs[0] << " unpaired=" << s_pairs_unmap[0] << " paired_other_ctg=" << s_pairs_ormap[0] << std::endl;
	std::cerr << "slaveMisEvid=" << slaveMisEvid << " masterMisEvid=" << masterMisEvid << std::endl;*/

	for( int i=0; i < m_mates.size(); i++ )
	{
		uint64_t masterTotPairs = m_mates[i] + m_mates_unmap[i] + m_mates_ormap[i];
		uint64_t slaveTotPairs = s_mates[i] + s_mates_unmap[i] + s_mates_ormap[i];
		
		if(slaveTotPairs >= 5 || masterTotPairs >= 5)
		{
			double masterScore = (masterTotPairs > 0) ? 100 * m_mates[i] / double(masterTotPairs) : 0.0;
			double slaveScore = (slaveTotPairs > 0) ? 100 * s_mates[i] / double(slaveTotPairs) : 0.0;
			
			if( masterScore >= 95.0 && slaveScore <= 5.0 ) slaveMisEvid++;
			if( slaveScore >= 95.0 && masterScore <= 5.0 ) masterMisEvid++;
		}
	}
    
    /*std::cerr << "MISASSEMBLY IN (" << masterId << "," << slaveId << ") MATES-STATS" << std::endl;
	std::cerr << "Master: paired=" << m_mates[0] << " unpaired=" << m_mates_unmap[0] << " paired_other_ctg=" << m_mates_ormap[0] << std::endl;
    std::cerr << "Slave:  paired=" << s_mates[0] << " unpaired=" << s_mates_unmap[0] << " paired_other_ctg=" << s_mates_ormap[0] << std::endl;
	std::cerr << "slaveMisEvid=" << slaveMisEvid << " masterMisEvid=" << masterMisEvid << std::endl;*/

	if( masterMisEvid > slaveMisEvid ) return MIS_MASTER;
	if( slaveMisEvid > masterMisEvid ) return MIS_SLAVE;

	return UNKNOWN; //return UNKNOWN;
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



PairedContig& PctgBuilder::addFirstContigTo(PairedContig& pctg, const int32_t ctgId) const
{
	const Contig& ctg = this->loadMasterContig(ctgId);

	if(ctg.size() == 0) return pctg;

	uint64_t start = 0;
	uint64_t end = ctg.size()-1;

	pctg.addMasterCtgId( ctgId );

	int64_t idx = pctg.size();
	uint64_t newBases = end - start + 1;
	pctg.resize( pctg.size() + newBases );

	for( int64_t i = start; i <= end; i++ ) pctg[idx++] = ctg.at(i);

	std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();
	mergeList.push_back( CtgInPctgInfo( ctgId, start, end, false, true ) );

	return pctg;

	/*if( pctg.size() != 0 ) throw std::invalid_argument( "Paired contig is not empty" );
    //PairedContig out(pctg.getId());
    const Contig& ctg = this->loadMasterContig(ctgId);

	ContigInPctgInfo ctgInfo(ctgId, ctg.size(), 0);
	pctg.getMasterCtgMap()[ctgInfo.getId()] = ctgInfo;

    return this->extendPctgWithCtgFrom(pctg,ctg,ctgInfo,std::pair<uint64_t,uint64_t>(0,0),true);*/
}


PairedContig& PctgBuilder::addFirstSlaveContigTo(PairedContig& pctg, const int32_t ctgId) const
{
	if( pctg.size() != 0 ) throw std::invalid_argument( "Paired contig is not empty" );
													   //PairedContig out(pctg.getId());
	const Contig& ctg = this->loadSlaveContig(ctgId);

	ContigInPctgInfo ctgInfo(ctgId, ctg.size(), 0);
	pctg.getSlaveCtgMap()[ctgInfo.getId()] = ctgInfo;

	return this->extendPctgWithCtgFrom(pctg,ctg,ctgInfo,std::pair<uint64_t,uint64_t>(0,0),true);
}


PairedContig& PctgBuilder::addFirstBlockTo(PairedContig& pctg, const std::list<Block> &blocks_list) const
{
	const Block &firstBlock = blocks_list.front();

    int32_t masterCtgId = firstBlock.getMasterId();
    pctg = addFirstContigTo(pctg,masterCtgId);

    return this->extendByBlock( pctg, blocks_list );
}


PairedContig& PctgBuilder::extendByBlock(PairedContig& pctg, const std::list<Block> &blocks_list) const
{
	if( pctg.size() == 0 ) return this->addFirstBlockTo( pctg, blocks_list );

    IdType masterCtgId, slaveCtgId;

	const Block &block = blocks_list.front();
	const Frame &master_frame = block.getMasterFrame();
	const Frame &slave_frame = block.getSlaveFrame();

    masterCtgId = master_frame.getContigId();
    slaveCtgId = slave_frame.getContigId();

    if( !pctg.containsMasterCtg(masterCtgId) )
    {
        if( !pctg.containsSlaveCtg(slaveCtgId) )
        {
            throw std::logic_error("Paired contig cannot be extended by this block.");
        }
        else
        {
            // merge the master contig of the block with the paired contig
            return this->mergeContig(pctg,blocks_list,true);
        }
    }
    else
    {
        if( !pctg.containsSlaveCtg(slaveCtgId) )
        {
            // merge the slave contig of the block with the paired contig
            return this->mergeContig(pctg,blocks_list,false);
        }
    }

    return pctg;
}


PairedContig&
PctgBuilder::mergeContig(
	PairedContig &pctg,
	const std::list<Block> &blocks_list,
	bool mergeMasterCtg ) const
{
	const Block &firstBlock = blocks_list.front();
	const Block &lastBlock = blocks_list.back();

	const Frame &firstMasterFrame = firstBlock.getMasterFrame();
	const Frame &firstSlaveFrame = firstBlock.getSlaveFrame();
	const Frame &lastMasterFrame = lastBlock.getMasterFrame();
	const Frame &lastSlaveFrame = lastBlock.getSlaveFrame();

    int32_t firstCtgId,  // identifier of the contig to be merged
			secondCtgId; // identifier of the contig already inside the paired contig

    if( mergeMasterCtg )
    {
        firstCtgId = firstMasterFrame.getContigId();
        secondCtgId = firstSlaveFrame.getContigId();
    }
    else
    {
        firstCtgId = firstSlaveFrame.getContigId();
        secondCtgId = firstMasterFrame.getContigId();
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
	//this->findBestAlignment(pctg,pctg.getContigInfo(secondCtgId,!mergeMasterCtg),startPos,endPos,ctg,mergeMasterCtg,blocks_list)
	BestPctgCtgAlignment bestAlign;

	if( bestAlign.main_homology() < MIN_HOMOLOGY )
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
	if( bestAlign.main_homology() < MIN_HOMOLOGY )
    {
// 		std::pair<uint64_t,uint64_t> masterId = std::make_pair( firstMasterFrame.getAssemblyId(), firstMasterFrame.getContigId() );
// 		std::pair<uint64_t,uint64_t> slaveId = std::make_pair( firstSlaveFrame.getAssemblyId(), firstSlaveFrame.getContigId() );

// 		double masterZscore = _tbp->computeZScore( masterId, (mergeMasterCtg ? ctg_beg : beginInCtg), (mergeMasterCtg ? ctg_end : endInCtg), true );
// 		double slaveZscore = _tbp->computeZScore( slaveId, (!mergeMasterCtg ? ctg_beg : beginInCtg), (!mergeMasterCtg ? ctg_end : endInCtg), false );

		pthread_mutex_lock(&g_badAlignMutex);

		if( blocks_list.size() < 1000 ) ext_fpi_2[ blocks_list.size() ]++;

		std::cerr << "bad main alignment (" << firstBlock.getMasterFrame().getContigId() << "," << firstBlock.getSlaveFrame().getContigId() << ")\n";
		std::cerr << "\t=> Start/End = pctg(" << startPos << "," << endPos << ") " << (mergeMasterCtg ? "master(" : "slave(") << ctg_beg << "," << ctg_end << ")\t";
		std::cerr << "Homology = " << bestAlign.main_homology() << "\tBegin = " << bestAlign[0].begin_a() << "/" << bestAlign[0].begin_b() << " Rev = " << bestAlign.isCtgReversed() << std::endl;
		std::cerr << "\t=> Regions = ([" << (mergeMasterCtg ? ctg_beg : beginInCtg) << "," << (mergeMasterCtg ? ctg_end : endInCtg) << "] -- ["
		<< (!mergeMasterCtg ? ctg_beg : beginInCtg) << "," << (!mergeMasterCtg ? ctg_end : endInCtg) << "])" << std::endl;
//		std::cerr << "\t=> Z-Scores = (" << masterZscore << "," << slaveZscore << ")" << std::endl;
		std::cerr << "\t=> Blocks = " << blocks_list.size() << std::endl;
		std::cerr << std::endl;

		pthread_mutex_unlock(&g_badAlignMutex);

		// adding description of possible duplicate region

		uint64_t con_evid = 0, dis_evid = 0;
		for( std::list<Block>::const_iterator b = blocks_list.begin(); b != blocks_list.end(); b++ )
		{
			const Frame& mf = b->getMasterFrame();
			const Frame& sf = b->getSlaveFrame();

			if( mf.getStrand() != sf.getStrand() ) dis_evid += b->getReadsNumber(); // discordant frames
			else con_evid += b->getReadsNumber(); // concordant frames
		}

		int64_t overlap_length = mergeMasterCtg ? (ctg_end - ctg_beg) : (endInCtg - beginInCtg);

		if( (double(con_evid) / double(con_evid+dis_evid)) >= 0.5 ) // ctg probably fwd
		{
			// compute pctg tails (i1,i2) and ctg tails (j1,j2)
			uint64_t i1 = startPos;
			uint64_t i2 = pctg.size() - endPos - 1;

			uint64_t j1 = mergeMasterCtg ? ctg_beg : beginInCtg;
			uint64_t j2 = ctg.size() - (mergeMasterCtg ? ctg_end : endInCtg) - 1;

			overlap_length += std::min(i1,j1) + std::min(i2,j2);
		}
		else // ctg probably reversed
		{
			// compute pctg tails (i1,i2) and ctg tails (j1,j2)
			uint64_t i1 = startPos;
			uint64_t i2 = pctg.size() - endPos - 1;

			uint64_t j1 = ctg.size() - (mergeMasterCtg ? ctg_end : endInCtg) - 1;
			uint64_t j2 = mergeMasterCtg ? ctg_beg : beginInCtg;

			overlap_length += std::min(i1,j1) + std::min(i2,j2);
		}

		pctg.addDupRegion( overlap_length, firstBlock.getSlaveFrame().getContigId() );

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
    else if( bestAlign.left().homology() < MIN_HOMOLOGY_2 || bestAlign.right().homology() < MIN_HOMOLOGY_2 ) // bad tail alignments
	{
		pthread_mutex_lock(&g_badAlignMutex);
		if( blocks_list.size() < 1000 )  ext_fpi_2[ blocks_list.size() ]++;
		//pthread_mutex_unlock(&g_badAlignMutex);

		//std::cerr << "good blocks alignment:\t(0," << firstBlock.getMasterFrame().getContigId() << ")\t("
		//	<< firstBlock.getSlaveFrame().getAssemblyId() << "," << firstBlock.getSlaveFrame().getContigId() << ")"
		//	<< "\tidentity=" << bestAlign.main_homology() << std::endl;

		int32_t masterId = firstMasterFrame.getContigId();
		int32_t slaveId = firstSlaveFrame.getContigId();
		double masterZscore = _tbp->computeZScore( masterBam, masterId, (mergeMasterCtg ? ctg_beg : beginInCtg), (mergeMasterCtg ? ctg_end : endInCtg), true );
		double slaveZscore = _tbp->computeZScore( slaveBam, slaveId, (!mergeMasterCtg ? ctg_beg : beginInCtg), (!mergeMasterCtg ? ctg_end : endInCtg), false );

		int32_t id1 = (mergeMasterCtg) ? firstBlock.getSlaveFrame().getContigId() : firstBlock.getMasterFrame().getContigId();
		int32_t id2 = (mergeMasterCtg) ? firstBlock.getMasterFrame().getContigId() : firstBlock.getSlaveFrame().getContigId();

		std::cerr << "bad tail alignment:\t" << (mergeMasterCtg ? "slave(" : "master(") << id1 << "," << beginInCtg << "," << endInCtg << ")"
			<< (!mergeMasterCtg ? "\tslave(" : "\tmaster(") << id2 << "," << ctg_beg << "," << ctg_end << ")";

		if( bestAlign.left().homology() < MIN_HOMOLOGY_2 ) std::cerr << "\tLEFT(" << bestAlign.left().homology() << ")";
		if( bestAlign.right().homology() < MIN_HOMOLOGY_2 ) std::cerr << "\tRIGHT(" << bestAlign.right().homology() << ")";
		std::cerr << std::endl;

		std::cerr << "z-scores = " << (mergeMasterCtg ? slaveZscore : masterZscore) << " / " << (!mergeMasterCtg ? slaveZscore : masterZscore)
			<< "\n" << std::endl;

		pthread_mutex_unlock(&g_badAlignMutex);

		// adding description of possible duplicate region

		uint64_t con_evid = 0, dis_evid = 0;
		for( std::list<Block>::const_iterator b = blocks_list.begin(); b != blocks_list.end(); b++ )
		{
			const Frame& mf = b->getMasterFrame();
			const Frame& sf = b->getSlaveFrame();

			if( mf.getStrand() != sf.getStrand() ) dis_evid += b->getReadsNumber(); // discordant frames
			else con_evid += b->getReadsNumber(); // concordant frames
		}

		int64_t overlap_length = mergeMasterCtg ? (ctg_end - ctg_beg) : (endInCtg - beginInCtg);

		if( (double(con_evid) / double(con_evid+dis_evid)) >= 0.5 ) // ctg probably fwd
		{
			// compute pctg tails (i1,i2) and ctg tails (j1,j2)
			uint64_t i1 = startPos;
			uint64_t i2 = pctg.size() - endPos - 1;

			uint64_t j1 = mergeMasterCtg ? ctg_beg : beginInCtg;
			uint64_t j2 = ctg.size() - (mergeMasterCtg ? ctg_end : endInCtg) - 1;

			overlap_length += std::min(i1,j1) + std::min(i2,j2);
		}
		else // ctg probably reversed
		{
			// compute pctg tails (i1,i2) and ctg tails (j1,j2)
			uint64_t i1 = startPos;
			uint64_t i2 = pctg.size() - endPos - 1;

			uint64_t j1 = ctg.size() - (mergeMasterCtg ? ctg_end : endInCtg) - 1;
			uint64_t j2 = mergeMasterCtg ? ctg_beg : beginInCtg;

			overlap_length += std::min(i1,j1) + std::min(i2,j2);
		}

		pctg.addDupRegion( overlap_length, firstBlock.getSlaveFrame().getContigId() );

		return pctg;
	}

	pthread_mutex_lock(&g_badAlignMutex);
	if( blocks_list.size() < 1000 )  ext_fpi[ blocks_list.size() ]++;
	pthread_mutex_unlock(&g_badAlignMutex);

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
		const int32_t ctgInPctgId,
        const int32_t ctgId,
        const BestPctgCtgAlignment& alignment,
        bool mergeMaster) const
{
	// following: previous code - DON'T DELETE
	/*if( mergeMaster ) return this->mergeMasterCtgInPos(pctg,ctg,ctgId,bestAlign);
	return this->mergeSlaveCtgInPos(pctg,ctg,ctgId,bestAlign);*/

	// new way of merging
	uint64_t pctg_shift = 0;
	std::pair<uint64_t,uint64_t> start_align, end_align;
	ContigInPctgInfo ctgInfo( ctgId, alignment );

	ContigInPctgInfo& ctgInPctgInfo = mergeMaster ? pctg.getSlaveCtgMap()[ctgInPctgId] : pctg.getMasterCtgMap()[ctgInPctgId];

	first_match_pos( alignment[0], start_align );
	last_match_pos( alignment[alignment.size()-1], end_align );

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
	uint64_t pctgScore = 0, ctgScore = 0;
	for( uint64_t i = start_align.first; i <= end_align.first; i++ ) if( char(pctg[i]) == 'N' ) pctgScore++;
	for( uint64_t i = start_align.second; i <= end_align.second; i++ ) if( char(ctg[i]) == 'N' ) ctgScore++;

	/*double pctgZscore = (this->_tbp)->computeZScore( ctgInPctgId, start_align.first, end_align.first, !mergeMaster );
	double ctgZscore = (this->_tbp)->computeZScore( ctgId, start_align.second, end_align.second, mergeMaster );

	if(pctgZscore < 0) pctgZscore = -pctgZscore;
	if(ctgZscore < 0) ctgZscore = -ctgZscore;*/

	if( ctgScore < pctgScore ) //if( ctg_n <= pctg_n )
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


void PctgBuilder::findBestAlignment(
        BestCtgAlignment &bestAlign,
        Contig &masterCtg,
		uint64_t masterStart,
		uint64_t masterEnd,
		Contig &slaveCtg,
        uint64_t slaveStart,
        uint64_t slaveEnd,
        const std::list<Block> &blocks_list ) const
{
    uint64_t con_evid = 0, dis_evid = 0;
	uint64_t mf_len = 0, sf_len = 0; // sum of master/slave frames lengths
	uint32_t blocks_num = blocks_list.size();

	int32_t min_frame_len = 100;

	int32_t masterId = -1, slaveId = -1;

	// compute the probability of slaveCtg to be reverse complemented respect to masterCtg
	for( std::list<Block>::const_iterator b = blocks_list.begin(); b != blocks_list.end(); b++ )
	{
		const Frame& mf = b->getMasterFrame();
		const Frame& sf = b->getSlaveFrame();

		int32_t min_len = std::min( mf.getLength(), sf.getLength() );

		if( b == blocks_list.begin() ) min_frame_len = min_len;
		else if( min_frame_len > min_len ) min_frame_len = min_len;

		masterId = mf.getContigId();
		slaveId = sf.getContigId();

		mf_len += mf.getLength();
		sf_len += sf.getLength();

		if( mf.getStrand() != sf.getStrand() ) dis_evid += b->getReadsNumber(); // discordant frames
			else con_evid += b->getReadsNumber(); // concordant frames
	}

	//uint64_t min_align_len = std::min(mf_len,sf_len);
	double con_prob = double(con_evid) / double(con_evid+dis_evid);

    // compute length thresholds
    size_t mt = 0.3 * masterCtg.size();
	size_t st = 0.3 * slaveCtg.size();

	int32_t align_threshold = 0.7 * min_frame_len;
	int32_t threshold = std::min( size_t(200), std::min(mt,st) );

	BandedSmithWaterman aligner; //BandedSmithWaterman aligner( bsw_band );
	MyAlignment align, bad_align(0), good_align(100);

	bool good_align_found = false;
	bool isSlaveRev = false;

	std::vector< MyAlignment > aligns;
    uint64_t tempPos;

	// contigs more likely have the same orientation
	if( con_prob >= 0.5 )
	{
		// align each block
		this->alignBlocks( masterCtg, masterStart, slaveCtg, slaveStart, blocks_list, aligns );
		
		// DEBUG
		/*for( size_t i=0; i < aligns.size(); i++ )
		{
			if( aligns[i].homology() < MIN_HOMOLOGY )
			{
				std::cerr << "[debug] bad alignment master(" << masterId << "," << aligns[i].begin_a() << "); "
				<< "slave(" << slaveId << "," << aligns[i].begin_b() << "); hom=" << aligns[i].homology() << std::endl;
			}
		}*/

		if( this->is_good( aligns, align_threshold ) )
		{
			good_align_found = true;
			isSlaveRev = false;
		}
		else
		{
			// else, try reversing the contig
			reverse_complement(slaveCtg);

			// update slave start/end positions
            tempPos = slaveStart;
            slaveStart = slaveCtg.size() - slaveEnd - 1;
            slaveEnd = slaveCtg.size() - tempPos - 1;

			// new alignments (one for each block)
			this->alignBlocks( masterCtg, masterStart, slaveCtg, slaveStart, blocks_list, aligns );

			// if align is good return it
			if( this->is_good( aligns, align_threshold ) )
            {
                good_align_found = true;
                isSlaveRev = true;
            }
		}
	}

	// contigs more likely have opposite orientations
	if( con_prob < 0.5 )
	{
		reverse_complement(slaveCtg);

		// update start and end positions of the blocks
		tempPos = slaveStart;
        slaveStart = slaveCtg.size() - slaveEnd - 1;
        slaveEnd = slaveCtg.size() - tempPos - 1;

		// new alignments (one for each block)
		this->alignBlocks( masterCtg, masterStart, slaveCtg, slaveStart, blocks_list, aligns );
		
		// DEBUG
		/*for( size_t i=0; i < aligns.size(); i++ )
		{
			if( aligns[i].homology() < MIN_HOMOLOGY )
			{
				std::cerr << "[debug] bad alignment master(" << masterId << "," << aligns[i].begin_a() << "); "
				<< "slave(" << slaveId << "," << aligns[i].begin_b() << "); hom=" << aligns[i].homology() << std::endl;
			}
		}*/

		if( this->is_good( aligns, align_threshold ) )
		{
			good_align_found = true;
			isSlaveRev = true;
		}
		else
		{
            // restore original (unreversed) contig
			reverse_complement(slaveCtg);

			// update start and end positions of the blocks
			tempPos = slaveStart;
            slaveStart = slaveCtg.size() - slaveEnd - 1;
            slaveEnd = slaveCtg.size() - tempPos - 1;

			// new alignments (one for each block)
			this->alignBlocks( masterCtg, masterStart, slaveCtg, slaveStart, blocks_list, aligns );

			if( this->is_good( aligns, align_threshold ) )
            {
                good_align_found = true;
                isSlaveRev = false;
            }
		}
	}

	// if the alignments computed were all bad, return a bad alignment to interrupt the merging
	if( !good_align_found || aligns.size() != blocks_num || blocks_num == 0 ){ bestAlign = BestCtgAlignment(bad_align,isSlaveRev); return; }

	// retrieve the starting/ending points of the alignment
	std::pair<uint64_t,uint64_t> alignStart, alignEnd;
	first_match_pos( aligns[0], alignStart );
	last_match_pos( aligns[blocks_num-1], alignEnd );

	// compute masterCtg tails (i1,i2) and slaveCtg tails (j1,j2)
	uint64_t i1 = alignStart.first;
	uint64_t i2 = masterCtg.size() - alignEnd.first - 1;
	uint64_t j1 = alignStart.second;
	uint64_t j2 = slaveCtg.size() - alignEnd.second - 1;

	// if tails are shorter than the computed threshold, assume the contigs have to be merged
    if( std::min(i1,j1) < threshold && std::min(i2,j2) < threshold ){ bestAlign = BestCtgAlignment(aligns,isSlaveRev); return; }

	MyAlignment leftAlign(100),rightAlign(100);
	bool leftRev, rightRev;

    ABlast ablast;

    /* LEFT TAIL ALIGNMENT */

    if( std::min(i1,j1) >= threshold )
	{
		if( i1 < j1 ) // masterCtg left tail < slaveCtg left tail
		{
            std::list<uint32_t> hits = ablast.findHits( slaveCtg, 0, alignStart.second-1, masterCtg, 0, alignStart.first-1 );
            size_t hits_num = hits.size();

            if( hits_num > 0 )
            {
                leftAlign = aligner.find_alignment( slaveCtg, hits.back(), alignStart.second-1, masterCtg, 0, alignStart.first-1, false, true );
                leftRev = true;
            }
            else
            {
                leftAlign = aligner.find_alignment( slaveCtg, alignStart.second-alignStart.first, alignStart.second-1, masterCtg, 0, alignStart.first-1, false, true );
                leftRev = true;
            }
		}
		else	// masterCtg left tail >= slaveCtg left tail
		{
            std::list<uint32_t> hits = ablast.findHits( masterCtg, 0, alignStart.first-1, slaveCtg, 0, alignStart.second-1 );
            size_t hits_num = hits.size();

            if( hits_num > 0 )
            {
                leftAlign = aligner.find_alignment( masterCtg, hits.back(), alignStart.first-1, slaveCtg, 0, alignStart.second-1, false, true );
                leftRev = false;
            }
            else
            {
                leftAlign = aligner.find_alignment( masterCtg, alignStart.first-alignStart.second, alignStart.first-1, slaveCtg, 0, alignStart.second-1, false, true );
                leftRev = false;
            }
		}
	}

    /* RIGHT TAIL ALIGNMENT */

	if( std::min(i2,j2) >= threshold )
	{
		if( i2 < j2 ) // pctg right tail < ctg right tail
		{
			Contig rightTail = chop_begin( slaveCtg, alignEnd.second+1 );

            std::list<uint32_t> hits = ablast.findHits( rightTail, 0, rightTail.size()-1, masterCtg, alignEnd.first+1, masterCtg.size()-1 );
            size_t hits_num = hits.size();

            if( hits_num > 0 )
            {
                rightAlign = aligner.find_alignment( rightTail, hits.front(), rightTail.size()-1, masterCtg, alignEnd.first+1, masterCtg.size()-1, true, false );
                rightRev = true;
            }
            else
            {
                rightAlign = aligner.find_alignment( rightTail, 0, rightTail.size()-1, masterCtg, alignEnd.first+1, masterCtg.size()-1, true, false );
                rightRev = true;
            }
		}
		else	// pctg right tail >= ctg right tail
		{
			Contig rightTail = chop_begin( masterCtg, alignEnd.first+1 );

            std::list<uint32_t> hits = ablast.findHits( rightTail, 0, rightTail.size()-1, slaveCtg, alignEnd.second+1, slaveCtg.size()-1 );
            size_t hits_num = hits.size();

            if( hits_num > 0 )
            {
                rightAlign = aligner.find_alignment( rightTail, hits.front(), rightTail.size()-1, slaveCtg, alignEnd.second+1, slaveCtg.size()-1, true, false );
                rightRev = false;
            }
            else
            {
                rightAlign = aligner.find_alignment( rightTail, 0, rightTail.size()-1, slaveCtg, alignEnd.second+1, slaveCtg.size()-1, true, false );
                rightRev = false;
            }
		}
	}

	bestAlign = BestCtgAlignment( aligns, isSlaveRev, leftAlign, rightAlign, leftRev, rightRev );
}


void PctgBuilder::alignBlocks(
	const Contig &masterCtg,
	const uint64_t &masterStart,
	const Contig &slaveCtg,
	const uint64_t &slaveStart,
	const std::list<Block> &blocks_list,
	std::vector< MyAlignment > &alignments ) const
{
	// initialize output
	alignments.clear();

	BandedSmithWaterman aligner; //ABlast aligner; //(this->_maxAlignment, this->_maxPctgGap, this->_maxCtgGap);

	// first & last blocks references
	const Block &firstBlock = blocks_list.front();
	const Block &lastBlock = blocks_list.back();

	// first & last frames (of pctg) references
	const Frame &masterFirstFrame = firstBlock.getMasterFrame();
	const Frame &masterLastFrame = lastBlock.getMasterFrame();

	int32_t masterId = firstBlock.getMasterId();
	int32_t slaveId = firstBlock.getSlaveId();

	int64_t masterStartAlign = masterStart;
	int64_t slaveStartAlign = slaveStart;
	uint32_t idx = 0;

	MyAlignment align;
	std::pair<uint64_t,uint64_t> last_match,first_match;

	Frame mf,prev_mf,sf,prev_sf;

	if( masterFirstFrame.getBegin() <= masterLastFrame.getBegin() ) // lista da processare in ordine
	{
		for( std::list<Block>::const_iterator b = blocks_list.begin(); b != blocks_list.end(); b++ )
		{
			mf = b->getMasterFrame();
			sf = b->getSlaveFrame();

			int32_t mlen = mf.getLength(); // lunghezza frame corrente su masterCtg
			int32_t slen = sf.getLength(); // lunghezza frame corrente su slaveCtg

			if( idx > 0 )
			{
				int32_t mgap = prev_mf.getBegin() <= mf.getBegin() ? (mf.getBegin() - prev_mf.getEnd() - 1) : (prev_mf.getBegin() - mf.getEnd() - 1);
				int32_t sgap = prev_sf.getBegin() <= sf.getBegin() ? (sf.getBegin() - prev_sf.getEnd() - 1) : (prev_sf.getBegin() - sf.getEnd() - 1);

				masterStartAlign = last_match.first + mgap; if( masterStartAlign < 0 ) masterStartAlign = 0;
				slaveStartAlign = last_match.second + sgap; if( slaveStartAlign < 0 ) slaveStartAlign = 0;
			}

			align = aligner.find_alignment( masterCtg, masterStartAlign, masterStartAlign+mlen-1, slaveCtg, slaveStartAlign, slaveStartAlign+slen-1 );
			alignments.push_back(align);

			last_match_pos( align, last_match );

			prev_mf = mf;
			prev_sf = sf;
			idx++;
		}
	}
	else // lista da processare in ordine inverso
	{
		for( std::list<Block>::const_reverse_iterator b = blocks_list.rbegin(); b != blocks_list.rend(); b++ )
		{
			mf = b->getMasterFrame();
			sf = b->getSlaveFrame();

			int32_t mlen = mf.getLength(); // lunghezza frame corrente su masterCtg
			int32_t slen = sf.getLength(); // lunghezza frame corrente su slaveCtg

			if( idx > 0 )
			{
				int32_t mgap = prev_mf.getBegin() <= mf.getBegin() ? (mf.getBegin() - prev_mf.getEnd() - 1) : (prev_mf.getBegin() - mf.getEnd() - 1);
				int32_t sgap = prev_sf.getBegin() <= sf.getBegin() ? (sf.getBegin() - prev_sf.getEnd() - 1) : (prev_sf.getBegin() - sf.getEnd() - 1);

				masterStartAlign = last_match.first + mgap; if( masterStartAlign < 0 ) masterStartAlign = 0;
				slaveStartAlign = last_match.second + sgap; if( slaveStartAlign < 0 ) slaveStartAlign = 0;
			}

			align = aligner.find_alignment( masterCtg, masterStartAlign, masterStartAlign+mlen-1, slaveCtg, slaveStartAlign, slaveStartAlign+slen-1 );
			alignments.push_back(align);

			last_match_pos( align, last_match );

			prev_mf = mf;
			prev_sf = sf;
			idx++;
		}
	}
}


bool PctgBuilder::is_good( const std::vector<MyAlignment> &aligns, uint64_t min_align_len ) const
{
	uint64_t align_len = 0;

	for( size_t i=0; i < aligns.size(); i++ )
	{
		if( aligns[i].homology() < MIN_HOMOLOGY ) return false;
		align_len += aligns[i].length();
	}

	return (align_len >= min_align_len);

	return true;//(align.homology() >= MIN_HOMOLOGY && align.length() >= min_align_len); // && align.score() > 0
}


bool PctgBuilder::is_good( const MyAlignment &align, uint64_t min_align_len ) const
{
	return (align.homology() >= MIN_HOMOLOGY && align.length() >= min_align_len); // && align.score() > 0
}

