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

#include <iostream>
#include <stdexcept>
#include <list>
#include <vector>

#include "alignment/my_alignment.hpp"

#define NUCLEOTIDE_PER_LINE 60

MyAlignment::MyAlignment() :
              _begin_a(0),
              _begin_b(0),
              _a_size(0),
              _b_size(0),
              _sequence(0),
              _score(0),
              _homology(0)
{}

MyAlignment::MyAlignment( const MyAlignment& orig ) :
              _begin_a(orig._begin_a),
              _begin_b(orig._begin_b),
              _a_size(orig._a_size),
              _b_size(orig._b_size),
              _sequence(orig._sequence),
              _score(orig._score),
			  _homology(orig._homology)
{}

MyAlignment::MyAlignment( double homology ):
		_begin_a(0),
		_begin_b(0),
		_a_size(0),
		_b_size(0),
		_sequence(0),
		_score(0),
		_homology(homology)
{}

MyAlignment::MyAlignment( size_type begin_a, size_type begin_b, size_type a_size, size_type b_size ) :
        _begin_a(begin_a), _begin_b(begin_b), _a_size(a_size), _b_size(b_size), _sequence(0), _homology(0)
{}

MyAlignment::MyAlignment(
	size_type begin_a,
	size_type begin_b,
	size_type a_size,
	size_type b_size,
	ScoreType score,
	double homology,
	const std::list<AlignmentAlphabet> edit_string ) :
		_begin_a(begin_a),
		_begin_b(begin_b),
		_a_size(a_size),
		_b_size(b_size),
		_score(score),
		_homology(homology)
{
    _sequence.reserve( edit_string.size() );

    std::list<AlignmentAlphabet>::const_iterator i;
    for( i = edit_string.begin(); i != edit_string.end(); i++ )
        _sequence.push_back( *i );
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

MyAlignment::size_type
MyAlignment::a_size() const
{
    return this->_a_size;
}

MyAlignment::size_type
MyAlignment::b_size() const
{
    return this->_b_size;
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
	return this->_homology;
	/*if( this->length() <= 0 ) return 0;

	uint64_t num_of_matches = 0;

    for( uint64_t i = 0; i < _sequence.size(); i++ )
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

    return ((RealType)num_of_matches) * 100 / ((RealType) this->length());*/
}

void
MyAlignment::set_homology( RealType homology )
{
    this->_homology = homology;
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
gaps_before_last_match( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &gaps )
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
    return this->a_position_in_b() + this->_a_size;
}

MyAlignment::int_type
MyAlignment::end_b_in_a() const
{
    return this->b_position_in_a() + this->_b_size;
}


void
MyAlignment::set_begin_a( size_type begin_a )
{
	this->_begin_a = begin_a;
}


void
MyAlignment::set_begin_b( size_type begin_b )
{
	this->_begin_b = begin_b;
}


const MyAlignment &
MyAlignment::operator=(const MyAlignment& orig)
{
    this->_begin_a = orig._begin_a;
    this->_begin_b = orig._begin_b;
    this->_a_size = orig._a_size;
    this->_b_size = orig._b_size;
    this->_score = orig._score;
    this->_sequence = orig._sequence;
	this->_homology = orig._homology;

    return *this;
}


void printAlignment( std::ostream& os, const Contig& a, const Contig& b, const MyAlignment& aln )
{
    os << "Contig A: "
       << "Name=\"" << a.name() << "\""
       << "\tLength=" << a.size() << std::endl;

    os << "Contig B: "
       << "Name=\"" << b.name() << "\""
       << "\tLength=" << b.size() << std::endl;

    os << std::endl;

    os << "Alignment: "
       //<< "Score=" << a.score()
       << "\tHomology=" << aln.homology()
       << "\tLength=" << aln.length()
       << "\tScore=" << aln.score()
       << std::endl;

    os << std::endl;

    unsigned int gap_num = 0, mismatch_num = 0;

    MyAlignment::size_type a_pos = aln.begin_a();
    MyAlignment::size_type b_pos = aln.begin_b();

    MyAlignment::SeqType::const_iterator i, begin_i;

    i = aln.sequence().begin();

    while( i != aln.sequence().end() )
    {
        if( i != aln.sequence().begin() ) os << std::endl;

        begin_i = i;
        unsigned int count = 0;

        os << a_pos << ":\t";
        while( count < NUCLEOTIDE_PER_LINE && i != aln.sequence().end() )
        {
            switch(*i)
            {
                case GAP_A:
                    os << "-";
                    break;
                case GAP_B:
                case MATCH:
                case MISMATCH:
                    os << a.at(a_pos);
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

        while( count < NUCLEOTIDE_PER_LINE && i != aln.sequence().end() )
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

        while( count < NUCLEOTIDE_PER_LINE && i != aln.sequence().end() )
        {
            switch(*i)
            {
                case GAP_B:
                    os << "-";
                    break;
                case GAP_A:
                case MATCH:
                case MISMATCH:
                    os << b.at(b_pos);
                    b_pos++;
                    break;
            }

            count++;
            i++;
        }

        os << std::endl;
    }

    std::cout << "\n\nGaps=" << gap_num << "\tMismatches=" << mismatch_num << "\n" << std::endl;
}