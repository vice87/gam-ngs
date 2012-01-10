#include <iostream>
#include <stdexcept>
#include <list>
#include <vector>

#include "alignment/my_alignment.hpp"

#define NUCLEOTIDE_PER_LINE 65 

MyAlignment::MyAlignment() : 
              _a(0),
              _b(0),
              _begin_a(0),
              _begin_b(0), 
              _sequence(0), 
              _score(0) 
{}

MyAlignment::MyAlignment( const MyAlignment& orig ) : 
              _a(orig._a),
              _b(orig._b),
              _begin_a(orig._begin_a),
              _begin_b(orig._begin_b), 
              _sequence(orig._sequence), 
              _score(orig._score) 
{}

MyAlignment::MyAlignment( const Contig &a, const Contig &b ) : 
        _a(a), _b(b), _begin_a(0), _begin_b(0), _sequence(0) 
{}

MyAlignment::MyAlignment( const Contig &a, size_type begin_a, const Contig &b, size_type begin_b ) :
        _a(a), _b(b), _begin_a(begin_a), _begin_b(begin_b), _sequence(0)
{}

MyAlignment::MyAlignment( const Contig &a, size_type begin_a, const Contig &b, size_type begin_b, ScoreType score, const std::list<AlignmentAlphabet> edit_string ) :
        _a(a), _b(b), _begin_a(begin_a), _begin_b(begin_b), _score(score)
{
    _sequence.reserve( edit_string.size() );
    
    std::list<AlignmentAlphabet>::const_iterator i;
    for( i = edit_string.begin(); i != edit_string.end(); i++ )
        _sequence.push_back( *i );
}

const Contig&
MyAlignment::a() const
{
    return this->_a;
}

const Contig&
MyAlignment::b() const
{
    return this->_b;
}

MyAlignment::size_type
MyAlignment::begin_a() const 
{
    return this->_begin_a; 
}

MyAlignment::size_type
MyAlignment::begin_b() const 
{
    return this->_begin_b; 
}

const MyAlignment::SeqType&
MyAlignment::sequence() const
{
    return this->_sequence;
}

ScoreType
MyAlignment::score() const
{
    return this->_score;
}

MyAlignment::size_type
MyAlignment::length() const
{
    return (this->_sequence).size();
}

RealType
MyAlignment::homology() const
{
    int_type num_of_matches = 0;
    
    for( MyAlignment::size_type i = 0; i < _sequence.size(); i++ ) 
    {
        switch(_sequence.at(i)) 
        {
            case MATCH:
                num_of_matches++;
                break;
            default:
                break;
        }
    }
    
    if( this->length() <= 0 ) return 0;
    
    return ((RealType)num_of_matches) * 100 / ((RealType) this->length());
}


bool
first_match_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos )
{
    pos.first = A.begin_a();
    pos.second = A.begin_b();
    
    for( MyAlignment::size_type i=0; i < A.sequence().size(); i++ )
    {
        switch( A.sequence().at(i) )
        {
            case MATCH:
                return true;
            case GAP_A:
                pos.second++;
                break;
            case GAP_B:
                pos.first++;
                break;
            case MISMATCH:
                pos.first++;
                pos.second++;
                break;
        }
    }
    
    return false;
}


bool
last_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos ) 
{
  pos.first = A.begin_a();
  pos.second = A.begin_b();
  bool at_least_one = false;

  for( MyAlignment::size_type i = 0; i < A.sequence().size(); i++ ) 
  {
      switch( A.sequence().at(i) )
      {
          case MATCH:
              at_least_one = true;
              pos.first++;
              pos.second++;
              break;
          case GAP_A:
              pos.second++;
              break;
          case GAP_B:
              pos.first++;
              break;
          case MISMATCH:
              pos.first++;
              pos.second++;
              break;
      }
  }
  
  return at_least_one;
}

bool
last_match_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos) 
{
    pos.first = A.begin_a();
    pos.second = A.begin_b();
    
    MyAlignment::size_type a_pos(A.begin_a()), b_pos(A.begin_b());
    bool at_least_one=false;
    
    for( MyAlignment::size_type i=0; i<A.sequence().size(); i++ ) 
    {
        switch(A.sequence().at(i)) 
        {
            case MATCH:
                at_least_one = true;
                pos.first = a_pos;
                pos.second = b_pos;
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
    
    return at_least_one;
}


bool
gaps_before_last_match( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> gaps ) 
{
    MyAlignment::size_type gaps_a_last(0), gaps_b_last(0);
    MyAlignment::size_type gaps_a(0), gaps_b(0);
    bool at_least_one = false;
    
    for( MyAlignment::size_type i=0; i< A.sequence().size(); i++ ) 
    {
        switch(A.sequence().at(i))
        {
            case MATCH:
                at_least_one = true;
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
    
    gaps.first = gaps_a_last;
    gaps.second = gaps_b_last;
    
    return at_least_one;
}


MyAlignment::int_type
MyAlignment::b_position_in_a() const 
{
    return (MyAlignment::int_type)this->begin_a() - (MyAlignment::int_type)this->begin_b();
}

MyAlignment::int_type
MyAlignment::a_position_in_b() const 
{
    return (MyAlignment::int_type)this->begin_b() - (MyAlignment::int_type)this->begin_a(); 
}


MyAlignment::int_type 
MyAlignment::end_a_in_b() const 
{
    return this->a_position_in_b()+this->_a.size(); 
}

MyAlignment::int_type 
MyAlignment::end_b_in_a() const 
{
    return this->b_position_in_a()+this->_b.size(); 
}


const MyAlignment &
MyAlignment::operator=(const MyAlignment& orig) 
{
    this->_a = orig._a;
    this->_b = orig._b;
    this->_begin_a = orig._begin_a;
    this->_begin_b = orig._begin_b;
    this->_score = orig._score;
    this->_sequence = orig._sequence;
    
    return *this;
}


std::ostream& operator<<( std::ostream& os, const MyAlignment& a ) 
{
    os << "Contig A: " 
       << "Name=\"" << a.a().name() << "\"" 
       << "\tLength=" << a.a().size() << std::endl;
    
    os << "Contig B: "
       << "Name=\"" << a.b().name() << "\"" 
       << "\tLength=" << a.b().size() << std::endl;
    
    os << std::endl;
    
    os << "Alignment: " 
       //<< "Score=" << a.score() 
       << "\tHomology=" << a.homology() 
       << "\tLength=" << a.length() 
       << "\tScore=" << a.score()
       << std::endl;
    
    os << std::endl;
    
    unsigned int gap_num = 0, mismatch_num = 0;
    
    MyAlignment::size_type a_pos = a.begin_a();
    MyAlignment::size_type b_pos = a.begin_b();
    
    MyAlignment::SeqType::const_iterator i, begin_i;
    
    i = a.sequence().begin();
    
    while( i != a.sequence().end() )
    {
        if( i != a.sequence().begin() ) os << std::endl;
        
        begin_i = i;
        unsigned int count = 0;
        
        os << a_pos << ":\t";
        while( count < NUCLEOTIDE_PER_LINE && i != a.sequence().end() )
        {
            switch(*i)
            {
                case GAP_A:
                    os << "-";
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
        
        i = begin_i;
        count = 0;
        
        while( count < NUCLEOTIDE_PER_LINE && i != a.sequence().end() )
        {
            switch(*i)
            {
                case GAP_A:
                case GAP_B:
                    os << " ";
                    gap_num++;
                    break;
                case MISMATCH:
                    os << " ";
                    mismatch_num++;
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
        
        i = begin_i;
        count = 0;
        
        while( count < NUCLEOTIDE_PER_LINE && i != a.sequence().end() )
        {
            switch(*i)
            {
                case GAP_B:
                    os << "-";
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
    
    std::cout << "\n\nGaps=" << gap_num << "\tMismatches=" << mismatch_num << "\n" << std::endl;
    
    return os;
}