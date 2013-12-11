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