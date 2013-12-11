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

#ifndef _MY_ALIGNMENT_
#define _MY_ALIGNMENT_

/*! \file
 *  \brief Definition of Alignment class.
 *  \details This file contains the defintion of Alignment.
 */

#include <list>
#include <vector>
#include <iostream>

#include "types.hpp"
#include "assembly/contig.hpp"

#define GAP_SCORE -8
#define MATCH_SCORE 5
#define MISMATCH_SCORE -4
#define GAP_EXT_SCORE -1

#ifndef FULL_ALIGNMENT_SIZE_PERCENTAGE
#define FULL_ALIGNMENT_SIZE_PERCENTAGE 95.00
#endif

typedef int64_t ScoreType;

typedef enum {
  GAP_A,
  GAP_B,
  MATCH,
  MISMATCH
} __attribute__((packed)) AlignmentAlphabet;


class MyAlignment
{
    friend void printAlignment( std::ostream& os, const Contig& a, const Contig& b, const MyAlignment& aln );

public:
    typedef int64_t int_type;
    typedef uint64_t size_type;
    typedef std::vector<AlignmentAlphabet> SeqType;

private:
    size_type _begin_a;
    size_type _begin_b;
    size_type _a_size;
    size_type _b_size;
    SeqType _sequence;
    ScoreType _score;
	double _homology;

public:

    MyAlignment();
    MyAlignment( const MyAlignment& orig );
	MyAlignment( double homology );
    MyAlignment( size_type begin_a, size_type begin_b, size_type a_size, size_type b_size );

	MyAlignment(
		size_type begin_a,
		size_type begin_b,
		size_type a_size,
		size_type b_size,
		ScoreType score,
		double homology,
		const std::list<AlignmentAlphabet> edit_string
	);

    size_type begin_a() const;
    size_type begin_b() const;

    size_type a_size() const;
    size_type b_size() const;

    const SeqType& sequence() const;

    size_type length() const;

    ScoreType score() const;

    RealType homology() const;

	int_type b_position_in_a() const;
    int_type a_position_in_b() const;

    int_type end_a_in_b() const;
    int_type end_b_in_a() const;

	void set_begin_a( size_type begin_a );
	void set_begin_b( size_type begin_b );

    void set_homology( RealType homology );

    const MyAlignment& operator=(const MyAlignment& orig);
};

bool first_match_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos );
bool last_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos );
bool last_match_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos );
bool gaps_before_last_match( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &gaps );

#endif // _MY_ALIGNMENT_
