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

#include <iostream>
#include <exception>
#include <boost/detail/container_fwd.hpp>
#include <sys/stat.h>

#include "pctg/PctgBuilder.hpp"
#include "pctg/MergeInCutTailFailed.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/ablast.hpp"
#include "alignment/banded_smith_waterman.hpp"

//std::ofstream g_badAlignStream;
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
}


PairedContig PctgBuilder::initByContig(const IdType& pctgId, const int32_t ctgId) const
{
    PairedContig pctg(pctgId);
    return this->addFirstContigTo(pctg,ctgId);
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


//TODO: split into several functions
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
	boost::tie(vbegin,vend) = boost::vertices(graph);
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

bool PctgBuilder::getMergePaths(const CompactAssemblyGraph &graph, CompactAssemblyGraph::Vertex &v, std::vector<MergeBlock> &mbv, MergeBlockLists &merge_paths) const {
    typedef CompactAssemblyGraph::Vertex Vertex;
    typedef CompactAssemblyGraph::VertexIterator VertexIterator;
    typedef CompactAssemblyGraph::OutEdgeIterator OutEdgeIterator;

    uint64_t out_degree = boost::out_degree(v, graph);
    uint64_t in_degree = boost::in_degree(v, graph);

    if (out_degree >= 2) {
        std::cerr << "[error] Found vertex with output degree >= 2 in fork-solved graph (this should NOT happen!) ==> (" << mbv[v].m_id << "," << mbv[v].s_id << ")" << std::endl;
        return false;
    }

    if (in_degree >= 2) {
        std::cerr << "[error] Found vertex with input degree >= 2 in fork-solved graph (this should NOT happen!) ==> (" << mbv[v].m_id << "," << mbv[v].s_id << ")" << std::endl;
        return false;
    }

    if (mbv[v].valid) merge_paths.front().push_back(mbv[v]);

    if (out_degree == 0) return true;

    if (out_degree == 1) {
        OutEdgeIterator ebegin, eend;
        boost::tie(ebegin, eend) = boost::out_edges(v, graph);

        Vertex v_cur = boost::source(*ebegin, graph);
        Vertex v_nxt = boost::target(*ebegin, graph);

        EdgeProperty edge_prop = boost::get(boost::edge_kind_t(), graph, *ebegin);
        double weight = edge_prop.weight;
        bool has_min_cov = edge_prop.min_cov;
        bool safe_edge = (weight >= 0 && weight < 0.3) || (weight < 0 && has_min_cov);

        //TODO: take care also of master edges, and edges with negative weight (i.e. not enough evidences)
        if (mbv[v_nxt].valid && edge_prop.kind == SLAVE_EDGE && safe_edge) {
            const std::list<Block> &cur_blocks = graph.getBlocks(v_cur);
            const std::list<Block> &nxt_blocks = graph.getBlocks(v_nxt);

            int32_t curSlaveStart = std::min(cur_blocks.front().getSlaveFrame().getBegin(), cur_blocks.back().getSlaveFrame().getBegin());
            int32_t nxtSlaveStart = std::min(nxt_blocks.front().getSlaveFrame().getBegin(), nxt_blocks.back().getSlaveFrame().getBegin());

            bool curFirstInSlave = curSlaveStart <= nxtSlaveStart;

            std::cerr << "[debug] Found (linear) SLAVE mis-assembly in ctg " << mbv[v_cur].s_id << std::endl;

            if (curFirstInSlave) {
                mbv[v_cur].s_rtail = false;
                mbv[v_nxt].s_ltail = false;
            } else {
                mbv[v_cur].s_ltail = false;
                mbv[v_nxt].s_rtail = false;
            }

            merge_paths.push_front(std::list<MergeBlock>());
            return this->getMergePaths(graph, v_nxt, mbv, merge_paths);
        }

        return this->getMergePaths(graph, v_nxt, mbv, merge_paths);
    }

    return false;
}


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

