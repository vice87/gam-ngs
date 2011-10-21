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
