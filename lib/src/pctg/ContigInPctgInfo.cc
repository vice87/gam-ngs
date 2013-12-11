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

#include "pctg/ContigInPctgInfo.hpp"
#include "types.hpp"
#include "pctg/BestPctgCtgAlignment.hpp"


ContigInPctgInfo::ContigInPctgInfo(): _reversed(false) {}

ContigInPctgInfo::ContigInPctgInfo(const int32_t ctgId, const UIntType& size, const UIntType& position):
	_reversed(false), _ctgId(ctgId), _size(size), _position(position), _left_cut(0), _right_cut(0) {}

ContigInPctgInfo::ContigInPctgInfo(const ContigInPctgInfo &orig):
	_ctgId(orig._ctgId), _size(orig._size),
	_position(orig._position), _reversed(orig._reversed),
	_left_cut(0), _right_cut(0) {}

ContigInPctgInfo::ContigInPctgInfo(const int32_t ctgId, const BestPctgCtgAlignment& bestAlign):
	_ctgId(ctgId), _left_cut(0), _right_cut(0)
{
    this->_size = bestAlign[0].b_size();
    this->_position = bestAlign[0].b_position_in_a();
    this->_reversed = bestAlign.isCtgReversed();
}


int32_t ContigInPctgInfo::getId() const
{
    return this->_ctgId;
}

uint64_t ContigInPctgInfo::getSize() const
{
	return this->_size;
}


int64_t ContigInPctgInfo::getFirstNucleotidePos() const
{
    return this->_position;
}

int64_t ContigInPctgInfo::getLastNucleotidePos() const
{
    return ((this->_position + this->_size > 0) ? (this->_position + this->_size - 1) : 0);
    //return this->_position + this->_size - this->_gaps - 1;
}

const bool& ContigInPctgInfo::isReversed() const
{
    return this->_reversed;
}

uint64_t ContigInPctgInfo::getLeftCut() const
{
	return this->_left_cut;
}

uint64_t ContigInPctgInfo::getRightCut() const
{
	return this->_right_cut;
}

void ContigInPctgInfo::setPosition(const UIntType& pos)
{
    this->_position = pos;
}

void ContigInPctgInfo::setLeftCut(const uint64_t& len)
{
	this->_left_cut = len;
}

void ContigInPctgInfo::setRightCut(const uint64_t& len)
{
	this->_right_cut = len;
}


const std::list< t_merge_gap >& ContigInPctgInfo::merge_gaps() const
{
	return this->_merge_gaps;
}

std::list< t_merge_gap >& ContigInPctgInfo::merge_gaps()
{
	return this->_merge_gaps;
}

void ContigInPctgInfo::addMergeGap( uint64_t start, uint64_t end, int64_t gap )
{
	t_merge_gap mg;

	mg.start = start;
	mg.end = end;
	mg.gap = gap;

	(this->_merge_gaps).push_back(mg);
}