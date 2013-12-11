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
