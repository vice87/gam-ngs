#include "pctg/CtgInPctgInfo.hpp"


CtgInPctgInfo::CtgInPctgInfo(): _reversed(false) {}

CtgInPctgInfo::CtgInPctgInfo(const int32_t ctgId, const int64_t start, const int64_t end, const bool reversed, const bool isMaster ) :
	_ctgId(ctgId), _start(start), _end(end), _reversed(reversed), _isMaster(isMaster)
{}

CtgInPctgInfo::CtgInPctgInfo(const CtgInPctgInfo &orig) :
	_ctgId(orig._ctgId), _start(orig._start), _end(orig._end), _reversed(orig._reversed), _isMaster(orig._isMaster)
{}


int32_t CtgInPctgInfo::getId() const
{
    return this->_ctgId;
}


int64_t CtgInPctgInfo::getStart() const
{
	return this->_start;
}


int64_t CtgInPctgInfo::getEnd() const
{
	return this->_end;
}


bool CtgInPctgInfo::isReversed() const
{
    return this->_reversed;
}

bool CtgInPctgInfo::isMaster() const
{
	return this->_isMaster;
}