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

#ifndef _ALIGNER_CODE_
#define _ALIGNER_CODE_

#include "alignment/aligner.hpp"

ScoreType
Aligner::match_score(const Contig& a, const Aligner::size_type& pos_a,
                     const Contig& b, const Aligner::size_type& pos_b) const {

  if (a.at(pos_a)==b.at(pos_b)) {
    return -2*(this->gap_score());
  }

  return 2*(this->_gap_score);
}

Aligner::Aligner(): _max_alignment(0), _gap_score(GAP_VALUE),
                                   _max_a_gaps(0),
                                   _max_b_gaps(0) {}

Aligner::Aligner(const ScoreType& gap_score):  _max_alignment(0),
                                   _gap_score(gap_score),
                                   _max_a_gaps(0),
                                   _max_b_gaps(0) {}

Aligner::Aligner(const Aligner::size_type& max_a_gaps,
                 const Aligner::size_type& max_b_gaps):
                                _max_alignment(0),
                                _gap_score(GAP_VALUE),
                                _max_a_gaps(max_a_gaps),
                                _max_b_gaps(max_b_gaps) {}

Aligner::Aligner(const ScoreType& gap_score,
                             const Aligner::size_type& max_a_gaps,
                             const Aligner::size_type& max_b_gaps):
                                    _max_alignment(0),
                                    _gap_score(gap_score),
                                    _max_a_gaps(max_a_gaps),
                                    _max_b_gaps(max_b_gaps) {}

Aligner::Aligner(const Aligner::size_type& max_alignment):
                                    _max_alignment(max_alignment),
                                    _gap_score(GAP_VALUE),
                                    _max_a_gaps(0),
                                    _max_b_gaps(0) {}

Aligner::Aligner(const Aligner::size_type& max_alignment,
                            const ScoreType& gap_score):
                                    _max_alignment(max_alignment),
                                    _gap_score(gap_score),
                                    _max_a_gaps(0),
                                    _max_b_gaps(0) {}

Aligner::Aligner(const Aligner::size_type& max_alignment,
                       const Aligner::size_type& max_a_gaps,
                       const Aligner::size_type& max_b_gaps):
                                    _max_alignment(max_alignment),
                                    _gap_score(GAP_VALUE),
                                    _max_a_gaps(max_a_gaps),
                                    _max_b_gaps(max_b_gaps) {}

Aligner::Aligner(const Aligner::size_type& max_alignment,
                       const ScoreType& gap_score,
                       const Aligner::size_type& max_a_gaps,
                       const Aligner::size_type& max_b_gaps):
                                    _max_alignment(max_alignment),
                                    _gap_score(gap_score),
                                    _max_a_gaps(max_a_gaps),
                                    _max_b_gaps(max_b_gaps) {}

ScoreType
Aligner::gap_score() const { return this->_gap_score; }

Aligner::size_type
Aligner::max_alignment() const { return this->_max_alignment; }

#endif // _ALIGNER_CODE_

