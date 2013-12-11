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

#ifndef _NUCLEOTIDE_CODE_
#define _NUCLEOTIDE_CODE_

#include <iostream>
#include <stdexcept>

#include "assembly/nucleotide.hpp"

Nucleotide::Nucleotide(): _base(N) {}

Nucleotide::Nucleotide(const Nucleotide& orig):
                               _base(orig._base) {}

Nucleotide::Nucleotide(const BaseType base_val):
                               _base(base_val) {}

Nucleotide::Nucleotide(const char base)
{
  switch (base) {
   case 'A':
   case 'a':
     this->_base=A;
     break;
   case 'T':
   case 't':
     this->_base=T;
     break;
   case 'C':
   case 'c':
     this->_base=C;
     break;
   case 'G':
   case 'g':
     this->_base=G;
     break;
   case 'N':
   case 'n':
     this->_base=N;
     break;
   default:
       this->_base=N;
       break;
     //throw std::invalid_argument("Unknown base value "+ std::string(&base) +" \n");
  }
}

const BaseType&
Nucleotide::base() const
{
  return _base;
}

const Nucleotide&
Nucleotide::operator=(const Nucleotide& orig)
{
  this->_base=orig._base;

  return *this;
}

const Nucleotide&
Nucleotide::operator=(const char orig)
{
  *this=Nucleotide(orig);

  return *this;
}

bool
Nucleotide::operator==(const Nucleotide& a)
{
  return this->_base==a._base;
}

bool
Nucleotide::operator!=(const Nucleotide& a)
{
  return this->_base!=a._base;
}

Nucleotide::operator char() const
{
  switch (this->_base) {
   case A:
     return 'A';
   case T:
     return 'T';
   case C:
     return 'C';
   case G:
     return 'G';
   case N:
   default:
     return 'N';
  }
}

Nucleotide
complement(const Nucleotide &orig)
{
  switch (orig._base) {
    case A:
      return Nucleotide(T);
    case T:
      return Nucleotide(A);
    case C:
      return Nucleotide(G);
    case G:
      return Nucleotide(C);
    case N:
    default:
      return Nucleotide(N);
  }
}

std::ostream& operator<<(std::ostream& os,
                         const Nucleotide& base)
{
  os << (char)base;

  return os;
}

#endif // _NUCLEOTIDE_
