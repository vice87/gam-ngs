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

#ifndef _ALIGNMENT_
#define _ALIGNMENT_

/*! \file
 *  \brief Definition of Alignment class.
 *  \details This file contains the defintion of Alignment.
 */

#ifndef FULL_ALIGNMENT_SIZE_PERCENTAGE
#define FULL_ALIGNMENT_SIZE_PERCENTAGE 95.00
#endif

#include "assembly/contig.hpp"

typedef long int ScoreType;

typedef enum {
  GAP_A,
  GAP_B,
  MATCH,
  MISMATCH
} __attribute__((packed)) AlignmentAlphabet;

class Alignment {
 public:
  typedef long int int_type;
  typedef unsigned long int size_type;
  typedef std::vector<AlignmentAlphabet> SeqType;
 private:
  Contig _a;
  Contig _b;
  size_type _begin_a;
  size_type _begin_b;
  bool _b_rev;
  SeqType _sequence;
  ScoreType _score;

 public:
  Alignment(const Contig &a, const Contig &b);

  Alignment(const Contig &a, const Contig &b,
                             const bool &b_rev);

  Alignment(const Contig &a, const Contig &b,
               const ScoreType score);

  /*Alignment(const Contig &a, const Contig &b,
            const bool &b_rev, const ScoreType score);*/

  Alignment(const Alignment& orig);

  Alignment(const Contig &a, const Contig &b,
            size_type begin_a, size_type begin_b);

  Alignment(const Contig &a, const Contig &b,
            size_type begin_a, size_type begin_b,
                                  const bool &b_rev);

  const size_type &begin_a() const;
  const size_type &begin_b() const;
  const bool &b_is_reversed() const;
  const SeqType &sequence() const;
  const size_type length() const;
  const ScoreType &score() const;
  const Contig &a() const;
  const Contig &b() const;

  bool is_full() const;

  bool is_full(const double& size_percentage) const;

  int_type end_a_in_b() const;

  int_type end_b_in_a() const;

  int_type b_position_in_a() const;

  int_type a_position_in_b() const;

  const Alignment &operator=(const Alignment& orig);

  double homology() const;

  const Alignment& a_gaps(const ScoreType &gap_score);

  const Alignment& b_gaps(const ScoreType &gap_score);

  const Alignment& match(const ScoreType &match_score);

  const Alignment& mismatch(const ScoreType &match_score);
};

Alignment a_gaps(const Alignment& orig, const ScoreType &gap_score);

Alignment b_gaps(const Alignment& orig, const ScoreType &gap_score);

Alignment match(const Alignment& orig, const ScoreType &match_score);

Alignment mismatch(const Alignment& orig, const ScoreType &match_score);

std::pair<Alignment::size_type,Alignment::size_type>
first_match_pos_in(const Alignment& A);

std::pair<Alignment::size_type,Alignment::size_type>
gaps_before_last_match_in(const Alignment& A);

std::pair<Alignment::size_type,Alignment::size_type>
last_pos_in(const Alignment& A);

std::pair<Alignment::size_type,Alignment::size_type>
last_match_pos_in(const Alignment& A);

std::ostream& operator<<(std::ostream& os, const Alignment& a);

#endif // _ALIGNMENT_
