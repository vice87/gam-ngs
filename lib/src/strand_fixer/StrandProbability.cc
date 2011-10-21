
#include "strand_fixer/StrandProbability.hpp"

StrandProbability::StrandProbability() :
        _probability( RealType(5)/RealType(10) ) {}


StrandProbability::StrandProbability( const StrandProbability& orig ) :
        _probability(orig._probability) {}


StrandProbability::StrandProbability( const RealType& probability ) :
        _probability(probability)
{
    this->boundValue();
}


void StrandProbability::boundValue() 
{
    if( this->_probability > 1 )
    {
        this->_probability = 1;
        return;
    }
    else if( this->_probability < 0 )
    {
        this->_probability = 0;
        return;
    }
}


const StrandProbability& StrandProbability::operator =( const StrandProbability &orig )
{
    this->_probability = orig._probability;
    return *this;
}


const StrandProbability& StrandProbability::operator =(const RealType& probability)
{
    this->_probability = probability;
    this->boundValue();
    return *this;
}


StrandProbability StrandProbability::operator *(const StrandProbability& sp) const
{
    StrandProbability output(*this);
    output._probability *= sp._probability;
    output.boundValue();
    
    return output;
}


StrandProbability StrandProbability::operator *(const RealType& prob) const
{
    StrandProbability output(*this);
    output._probability *= prob;
    output.boundValue();
    
    return output;
}


char StrandProbability::getStrand() const
{
    RealType half = RealType(5) / RealType(10);
    
    if( *this < half ) return '-';
    if( *this > half ) return '+';
    
    return '?';
}


StrandProbability::operator  RealType() const
{
    return this->_probability;
}