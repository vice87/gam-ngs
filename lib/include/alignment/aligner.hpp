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

#ifndef _ALIGNER_
#define _ALIGNER_

#include "alignment/alignment.hpp"

#define GAP_VALUE -1

class Aligner {
  public:
     typedef Alignment::int_type int_type;
     typedef Alignment::size_type size_type;

  protected:
     const size_type _max_alignment;
     const ScoreType _gap_score;

     size_type _max_a_gaps;
     size_type _max_b_gaps;

     ScoreType
     match_score(const Contig& a, const size_type& pos_a,
                 const Contig& b, const size_type& pos_b) const;
  public:

     Aligner();

     Aligner(const ScoreType& gap_score);

     Aligner(const size_type& max_a_gaps, const size_type& max_b_gaps);

     Aligner(const ScoreType& gap_score, const size_type& max_a_gaps,
                const size_type& max_b_gaps);

     Aligner(const size_type& max_alignment);

     Aligner(const size_type& max_alignment, const ScoreType& gap_score);

     Aligner(const size_type& max_alignment, const size_type& max_a_gaps,
                                      const size_type& max_b_gaps);

     Aligner(const size_type& max_alignment, const ScoreType& gap_score,
                const size_type& max_a_gaps, const size_type& max_b_gaps);

     virtual ~Aligner() {}

     ScoreType gap_score() const;

     size_type max_alignment() const;

     virtual Alignment
     find_alignment(const Contig& a, const size_type& begin_a,
                    const Contig& b, const size_type& begin_b) const = 0;
};

#endif // _ALIGNER_

