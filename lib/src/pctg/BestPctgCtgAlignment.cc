#include "pctg/BestPctgCtgAlignment.hpp"


BestPctgCtgAlignment::BestPctgCtgAlignment(const MyAlignment& a, const bool isCtgReverse):
        _a(a), _isCtgReverse(isCtgReverse) {}

BestPctgCtgAlignment::BestPctgCtgAlignment(const BestPctgCtgAlignment& orig):
        _a(orig._a), _isCtgReverse(orig._isCtgReverse) {}

const MyAlignment&
BestPctgCtgAlignment::getAlignment() const
{
    return this->_a;
}

const bool&
BestPctgCtgAlignment::isCtgReversed() const
{
    return this->_isCtgReverse;
}
		
const BestPctgCtgAlignment&
BestPctgCtgAlignment::operator =(const BestPctgCtgAlignment& orig)
{
    this->_a = (orig._a);
    this->_isCtgReverse = (orig._isCtgReverse);
    
    return *this;
}