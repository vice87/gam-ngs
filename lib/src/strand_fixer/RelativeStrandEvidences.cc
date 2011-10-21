#include "strand_fixer/RelativeStrandEvidences.hpp"


RelativeStrandEvidences::RelativeStrandEvidences() :
        _positive(0), _negative(0) {}


RelativeStrandEvidences::RelativeStrandEvidences(const RelativeStrandEvidences& orig):
        _positive(orig._positive), _negative(orig._negative) {}


const UIntType&
RelativeStrandEvidences::getPositiveEvidences() const
{
    return this->_positive;
}


const UIntType&
RelativeStrandEvidences::getNegativeEvidences() const
{
    return this->_negative;
}


UIntType
RelativeStrandEvidences::getEvidences() const
{
    return this->_positive + this->_negative;
}


void
RelativeStrandEvidences::addPositiveEvidences(const UIntType& positive)
{
    this->_positive += positive;
}


void
RelativeStrandEvidences::addNegativeEvidences(const UIntType& negative)
{
    this->_negative += negative;
}


const RelativeStrandEvidences& 
RelativeStrandEvidences::operator =(const RelativeStrandEvidences &orig)
{
    this->_positive = orig._positive;
    this->_negative = orig._negative;
    
    return *this;
}


RelativeStrandEvidences
RelativeStrandEvidences::operator +(const RelativeStrandEvidences &orig)
{
    RelativeStrandEvidences output;
    
    output._positive = this->_positive + orig._positive;
    output._negative = this->_negative + orig._negative;
    
    return output;
}


const RelativeStrandEvidences&
RelativeStrandEvidences::operator +=(const RelativeStrandEvidences &orig)
{
    this->_positive += orig._positive;
    this->_negative += orig._negative;
    
    return *this;
}


char 
RelativeStrandEvidences::getStrand() const
{
    if( this->_positive < this->_negative ) return '-';
    
    if( this->_positive > this->_negative ) return '+';
    
    return '?';
}
