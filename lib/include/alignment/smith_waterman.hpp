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

#ifndef _SMITH_WATERMAN_
#define _SMITH_WATERMAN_

#include "alignment/aligner.hpp"

class SmithWaterman: public Aligner {
   public:
     typedef Aligner::int_type int_type;
     typedef Aligner::size_type size_type;

   private:
     inline
     size_type get_b_pos(const size_type &y, const size_type &begin_b) const
     {
       return y+begin_b;
     }

     inline
     size_type get_a_pos(const size_type &y, const size_type &x,
                                          const size_type &begin_a) const
     {
       return x+y+begin_a-this->_max_b_gaps;
     }

     inline
     size_type get_y_pos(const size_type &b_pos,
                         const size_type &begin_b) const
     {
       return b_pos-begin_b;
     }

     inline
     size_type get_x_pos(const size_type &y, const size_type &a_pos,
                                          const size_type &begin_a) const
     {
       return a_pos+this->_max_b_gaps-begin_a-y;
     }

     void
     compute_alignment_matrix(const Contig& a, size_type& begin_a,
                              const Contig& b, size_type& begin_b,
                              ScoreType* a_matrix[],
                              const size_type& y_size,
                              const size_type& x_size) const;

     Alignment
     compute_alignment_from_matrix(const Contig& a, size_type& begin_a,
                              const Contig& b, size_type& begin_b,
                              ScoreType* a_matrix[],
                              const size_type& y_size,
                              const size_type& x_size,
                              const bool& b_rev) const;
  public:
     SmithWaterman();

     SmithWaterman(const size_type& max_a_gaps,
                   const size_type& max_b_gaps);

     SmithWaterman(const ScoreType& gap_score,
                                   const size_type& max_a_gaps,
                                   const size_type& max_b_gaps);

     SmithWaterman(const size_type& max_alignment);

     SmithWaterman(const size_type& max_alignment,
                                   const size_type& max_a_gaps,
                                   const size_type& max_b_gaps);

     SmithWaterman(const size_type& max_alignment,
                                   const ScoreType& gap_score,
                                   const size_type& max_a_gaps,
                                   const size_type& max_b_gaps);

     Alignment
     apply(const Contig& a, size_type begin_a,
            const Contig& b, size_type begin_b, const bool& b_rev) const;

     Alignment
     apply(const Contig& a, size_type begin_a,
            const Contig& b, size_type begin_b) const;

     inline
     Alignment
     find_alignment(const Contig& a, const size_type& begin_a,
                     const Contig& b, const size_type& begin_b) const
     {
       return apply(a,begin_a,b,begin_b,false);
     }
};

#endif // _SMITH_WATERMAN_

