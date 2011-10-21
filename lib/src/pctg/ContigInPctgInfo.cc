#include "pctg/ContigInPctgInfo.hpp"
#include "types.hpp"
#include "pctg/BestPctgCtgAlignment.hpp"


ContigInPctgInfo::ContigInPctgInfo(): _gaps(0), _reversed(false) {}

ContigInPctgInfo::ContigInPctgInfo(const IdType& ctgId, const UIntType& size, const UIntType& position):
        _gaps(0), _reversed(false), _ctgId(ctgId), _size(size), _position(position) {}

ContigInPctgInfo::ContigInPctgInfo(const ContigInPctgInfo &orig):
        _ctgId(orig._ctgId), _size(orig._size), _position(orig._position),
        _gaps(orig._gaps), _reversed(orig._reversed) {}

ContigInPctgInfo::ContigInPctgInfo(const IdType& ctgId, const BestPctgCtgAlignment& bestAlign):
        _ctgId(ctgId), _gaps(0)
{
    this->_size = bestAlign.getAlignment().b().size();
    this->_position = bestAlign.getAlignment().b_position_in_a();
    this->_reversed = bestAlign.isCtgReversed();
}


IdType ContigInPctgInfo::getId() const
{
    return this->_ctgId;
}

UIntType ContigInPctgInfo::getFirstNucleotidePos() const
{
    return this->_position;
}

UIntType ContigInPctgInfo::getLastNucleotidePos() const
{
    return this->_position + this->_size - this->_gaps - 1;
}

const bool& ContigInPctgInfo::isReversed() const
{
    return this->_reversed;
}

void ContigInPctgInfo::setGaps(const UIntType& gaps)
{
    this->_gaps = gaps;
}

void ContigInPctgInfo::setPosition(const UIntType& pos)
{
    this->_position = pos;
}