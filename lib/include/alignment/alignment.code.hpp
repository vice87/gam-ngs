#ifndef _ALIGNMENT_CODE_
#define _ALIGNMENT_CODE_

#define MIN_SCORE -10000000
#define NUCLEOTIDE_PER_LINE 65 

#include <iostream>
#include <stdexcept>

#include "alignment/alignment.hpp"

Alignment::Alignment(const Contig &a, const Contig &b): _a(a), _b(b),
         _begin_a(0), _begin_b(0), _b_rev(false), _sequence(0), _score(0) {}

Alignment::Alignment(const Contig &a, const Contig &b, const bool& b_rev): _a(a), _b(b),
         _begin_a(0), _begin_b(0), _b_rev(b_rev), _sequence(0), _score(0) {}

Alignment::Alignment(const Contig &a, const Contig &b, 
               const ScoreType score): _a(a), _b(b),
                               _begin_a(0), _begin_b(0),
                               _b_rev(false),
                               _sequence(0), _score(score) {}

/*Alignment::Alignment(const Contig &a, const Contig &b, const bool& b_rev, 
               const ScoreType score): _a(a), _b(b),
                               _begin_a(0), _begin_b(0),
                               _b_rev(b_rev),
                               _sequence(0), _score(score) {}*/

Alignment::Alignment(const Alignment& orig): 
              _a(orig._a),
              _b(orig._b),
              _begin_a(orig._begin_a),
              _begin_b(orig._begin_b), 
              _b_rev(orig._b_rev),
              _sequence(orig._sequence), 
              _score(orig._score) {}

Alignment::Alignment(const Contig &a, const Contig &b, 
            size_type begin_a, size_type begin_b): 
              _a(a), _b(b),
              _begin_a(begin_a), _begin_b(begin_b),
              _b_rev(false),
              _sequence(0), _score(0) {}
 
Alignment::Alignment(const Contig &a, const Contig &b, 
            size_type begin_a, size_type begin_b, 
                                const bool& b_rev): 
              _a(a), _b(b),
              _begin_a(begin_a), _begin_b(begin_b),
              _b_rev(b_rev),
              _sequence(0), _score(0) {}

const Alignment::size_type &
Alignment::begin_a() const 
{ 
  return this->_begin_a; 
}

const Alignment::size_type &
Alignment::begin_b() const 
{ 
  return this->_begin_b; 
}

const bool&
Alignment::b_is_reversed() const 
{ 
  return this->_b_rev; 
}

const Alignment::SeqType &
Alignment::sequence() const 
{ 
  return this->_sequence; 
}

const Alignment::size_type 
Alignment::length() const 
{ 
  return (this->_sequence).size(); 
}

const ScoreType &
Alignment::score() const 
{ 
  return this->_score; 
}

const Contig &
Alignment::a() const 
{ 
  return this->_a; 
}

const Contig &
Alignment::b() const 
{ 
  return this->_b; 
}

Alignment::int_type 
Alignment::end_a_in_b() const 
{ 
  return this->a_position_in_b()+this->_a.size(); 
}

Alignment::int_type 
Alignment::end_b_in_a() const 
{ 
  return this->b_position_in_a()+this->_b.size(); 
}

Alignment::int_type
Alignment::b_position_in_a() const 
{ 
  return (Alignment::int_type)this->begin_a()-
                 (Alignment::int_type)this->begin_b(); 
}

Alignment::int_type
Alignment::a_position_in_b() const 
{ 
  return (int)this->begin_b()-(int)this->begin_a(); 
}

const Alignment &
Alignment::operator=(const Alignment& orig) 
{
  this->_a=orig.a();
  this->_b=orig.b();
  this->_begin_a=orig.begin_a();
  this->_begin_b=orig.begin_b();
  this->_b_rev=orig.b_is_reversed();
  this->_sequence=orig.sequence();
  this->_score=orig.score();

  return *this;
}

bool
Alignment::is_full(const double& size_percentage) const
{
  try {
    std::pair<size_type,size_type> last_pos=
                     last_pos_in(*this);

    if (b_position_in_a()>=0)  {
      //if (last_pos.second<a().size()) {
        if (last_pos.second-_begin_b >= 10) {
          return true;
        } else {
          return false;
        }
      /*} else {
        if (((double)100)*(last_pos.second-_begin_b)>=
              size_percentage*last_pos.second) {
          return true;
        } else {
          return false;
        }
      }*/
    } else {
      //if (last_pos.first<b().size()) {
        if (last_pos.first-_begin_a >= 10) {
          return true;
        } else {
          return false;
        }
      /*} else {
        if (((double)100)*(last_pos.first-_begin_a)>=
               size_percentage*last_pos.first) {
          return true;
        } else {
          return false;
        }
      }*/
    }
  } catch (std::invalid_argument& e) {
    return false;
  }
}

bool
Alignment::is_full() const
{
  return this->is_full(FULL_ALIGNMENT_SIZE_PERCENTAGE);
}

double
Alignment::homology() const
{
  int_type num_of_matches=0;
  
  for (Alignment::size_type i=0; i<_sequence.size(); i++) {
    switch(_sequence.at(i)) {
      case MATCH:
        num_of_matches++;
        break;
      default:
        break;
    }
  }

  return ((double)num_of_matches)*100/((double) this->length());
}

const Alignment&
Alignment::a_gaps(const ScoreType &gap_score) 
{

  this->_sequence.push_back(GAP_A); 
  this->_score+=gap_score;

  return *this;
}

const Alignment&
Alignment::b_gaps(const ScoreType &gap_score) 
{

  this->_sequence.push_back(GAP_B); 
  this->_score+=gap_score;

  return *this;
}

const Alignment&
Alignment::match(const ScoreType &match_score) 
{

  this->_sequence.push_back(MATCH); 
  this->_score+=match_score;

  return *this;
}

const Alignment&
Alignment::mismatch(const ScoreType &match_score) 
{

  this->_sequence.push_back(MISMATCH); 
  this->_score+=match_score;

  return *this;
}

Alignment
a_gaps(const Alignment& orig, const ScoreType &gap_score) 
{
  Alignment new_a(orig);
  
  return new_a.a_gaps(gap_score);
}
  
Alignment
b_gaps(const Alignment& orig, const ScoreType &gap_score)
{
  Alignment new_a(orig);
  
  return new_a.b_gaps(gap_score);
}
  
Alignment
match(const Alignment& orig, const ScoreType &match_score) 
{
  Alignment new_a(orig);
  
  return new_a.match(match_score);
}

Alignment
mismatch(const Alignment& orig, const ScoreType &match_score) 
{
  Alignment new_a(orig);
  
  return new_a.mismatch(match_score);
}

std::pair<Alignment::size_type,Alignment::size_type>
first_match_pos_in(const Alignment& A) 
{
  Alignment::size_type a_pos(A.begin_a()),b_pos(A.begin_b());

  for (Alignment::size_type i=0; i<A.sequence().size(); i++) {
    switch(A.sequence().at(i)) {
      case MATCH:
        return 
           std::pair<Alignment::size_type,Alignment::size_type>(a_pos,b_pos);
      case GAP_A:
        b_pos++;
        break;
      case GAP_B:
        a_pos++;
        break;
      case MISMATCH:
        a_pos++;
        b_pos++;
        break;
    }
  }

  throw std::invalid_argument("No match in alignment\n");
}

std::pair<Alignment::size_type,Alignment::size_type>
gaps_before_last_match_in(const Alignment& A) 
{
  Alignment::size_type gaps_a_last(0), gaps_b_last(0);
  Alignment::size_type gaps_a(0), gaps_b(0);
  bool at_least_one=false;

  for (Alignment::size_type i=0; i<A.sequence().size(); i++) {
    switch(A.sequence().at(i)) {
      case MATCH:
        at_least_one=true;
        gaps_a_last=gaps_a;
        gaps_b_last=gaps_b;
        break;
      case GAP_A:
        gaps_a++;
        break;
      case GAP_B:
        gaps_b++;
        break;
      case MISMATCH:
        break;
    }
  }

  if (at_least_one) 
      return std::pair<Alignment::size_type,
                       Alignment::size_type>(gaps_a_last,gaps_b_last);

  throw std::invalid_argument("No match in alignment\n");
}

std::pair<Alignment::size_type,Alignment::size_type>
last_pos_in(const Alignment& A) 
{
  Alignment::size_type a_pos(A.begin_a()),b_pos(A.begin_b());
  bool at_least_one=false;

  for (Alignment::size_type i=0; i<A.sequence().size(); i++) {
    switch(A.sequence().at(i)) {
      case MATCH:
        at_least_one=true;
        a_pos++;
        b_pos++;
        break;
      case GAP_A:
        b_pos++;
        break;
      case GAP_B:
        a_pos++;
        break;
      case MISMATCH:
        a_pos++;
        b_pos++;
        break;
    }
  }
  
  if (at_least_one) 
    return std::pair<Alignment::size_type,Alignment::size_type>(a_pos,b_pos);

  throw std::invalid_argument("No match in alignment\n");
}

std::pair<Alignment::size_type,Alignment::size_type>
last_match_pos_in(const Alignment& A) 
{
  Alignment::size_type a_pos(A.begin_a()),b_pos(A.begin_b());
  Alignment::size_type a_last(a_pos), b_last(b_pos);
  bool at_least_one=false;

  for (Alignment::size_type i=0; i<A.sequence().size(); i++) {
    switch(A.sequence().at(i)) {
      case MATCH:
        at_least_one=true;
        a_last=a_pos;
        b_last=b_pos;
        a_pos++;
        b_pos++;
        break;
      case GAP_A:
        b_pos++;
        break;
      case GAP_B:
        a_pos++;
        break;
      case MISMATCH:
        a_pos++;
        b_pos++;
        break;
    }
  }
  
  if (at_least_one) 
    return std::pair<Alignment::size_type,Alignment::size_type>(a_last,b_last);

  throw std::invalid_argument("No match in alignment\n");
}

std::ostream& operator<<(std::ostream& os, const Alignment& a) 
{
  os << "Contig A: " 
     << "Name=\"" << a.a().name() << "\"" 
     << "\tLength=" << a.a().size() << "\tStrand=+"<< std::endl;
  
  os << "Contig B: " 
     << "Name=\"" << a.b().name() << "\"" 
     << "\tLength=" << a.b().size() << "\tStrand=";

  if (a.b_is_reversed()) {
    os << "-" << std::endl; 
  } else {
    os << "+" << std::endl;
  }

  os << std::endl;

  os << "Alignment: " 
     << "Score=" << a.score() 
     << "\tHomology=" << a.homology() 
     << "\tLength=" << a.length() << std::endl;

  os << std::endl;

  Alignment::size_type a_pos=a.begin_a();
  Alignment::size_type b_pos=a.begin_b();

  std::vector<AlignmentAlphabet>::const_iterator i,begin_i;

  i=a.sequence().begin();

  while (i!=a.sequence().end()) {

    if (i!=a.sequence().begin()) {
      os << std::endl;
    }

    begin_i=i;
    unsigned int count=0;

    os << a_pos << ":\t";
    while ((count < NUCLEOTIDE_PER_LINE)&&(i!=a.sequence().end())) {

      switch(*i){
       case GAP_A:
        os << "n";
        break;
       case GAP_B:
       case MATCH:
       case MISMATCH:
        os << a.a().at(a_pos);
        a_pos++;
        break;
      } 
      count++;
      i++;
    }
    os << std::endl;
    os << "\t";

    i=begin_i;
    count=0;
    while ((count < NUCLEOTIDE_PER_LINE)&&(i!=a.sequence().end())) {

      switch(*i){
       case GAP_A:
       case GAP_B:
        os << " ";
        break;
       case MISMATCH:
        os << "!";
        break;
       case MATCH:
        os << "|";
        break;
      } 
      count++;
      i++;
    }
    os << std::endl;

    os << b_pos << ":\t";
    i=begin_i;
    count=0;
    while ((count < NUCLEOTIDE_PER_LINE)&&(i!=a.sequence().end())) {

      switch(*i){
       case GAP_B:
        os << "n";
        break;
       case GAP_A:
       case MATCH:
       case MISMATCH:
        os << a.b().at(b_pos);
        b_pos++;
        break;
      } 
      count++;
      i++;
    }
    os << std::endl;
  }

  return os;
}

#endif // _ALIGNMENT_CODE_
