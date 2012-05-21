#include "pctg/ContigInPctgInfo.hpp"
#include "types.hpp"
#include "pctg/BestPctgCtgAlignment.hpp"


ContigInPctgInfo::ContigInPctgInfo(): _reversed(false) {}

ContigInPctgInfo::ContigInPctgInfo(const IdType& aId, const IdType& ctgId, const UIntType& size, const UIntType& position):
	_reversed(false), _aId(aId), _ctgId(ctgId), _size(size), _position(position), _left_cut(0), _right_cut(0) {}

ContigInPctgInfo::ContigInPctgInfo(const IdType& ctgId, const UIntType& size, const UIntType& position):
	_reversed(false), _ctgId(ctgId), _size(size), _position(position), _left_cut(0), _right_cut(0) {}

ContigInPctgInfo::ContigInPctgInfo(const ContigInPctgInfo &orig):
	_aId(orig._aId),_ctgId(orig._ctgId), _size(orig._size),
	_position(orig._position), _reversed(orig._reversed),
	_left_cut(0), _right_cut(0) {}

ContigInPctgInfo::ContigInPctgInfo(const IdType& aId, const IdType& ctgId, const BestPctgCtgAlignment& bestAlign):
	_aId(aId), _ctgId(ctgId), _left_cut(0), _right_cut(0)
{
	std::pair< uint64_t, uint64_t > start_align;
	first_match_pos( bestAlign[0], start_align );

    this->_size = bestAlign[0].b_size();
	this->_position = int64_t(start_align.first) - int64_t(start_align.second);
    this->_reversed = bestAlign.isCtgReversed();
}

ContigInPctgInfo::ContigInPctgInfo(const IdType& ctgId, const BestPctgCtgAlignment& bestAlign):
	_ctgId(ctgId), _left_cut(0), _right_cut(0)
{
    this->_size = bestAlign[0].b_size();
    this->_position = bestAlign[0].b_position_in_a();
    this->_reversed = bestAlign.isCtgReversed();
}


std::pair<IdType,IdType> ContigInPctgInfo::getId() const
{
    return std::make_pair(this->_aId,this->_ctgId);
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