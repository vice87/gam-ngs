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

#include "alignment/alignment.code.hpp"

class Alignment;

Alignment a_gaps(const Alignment& orig, const ScoreType &gap_score);

Alignment b_gaps(const Alignment& orig, const ScoreType &gap_score);

Alignment match(const Alignment& orig, const ScoreType &match_score);

Alignment mismatch(const Alignment& orig, const ScoreType &match_score);

std::pair<Alignment::size_type,Alignment::size_type>
first_match_pos_in(const Alignment& A);

std::pair<Alignment::size_type,Alignment::size_type>
gaps_before_last_match_in(const Alignment& A);

std::pair<Alignment::size_type,Alignment::size_type>
last_match_pos_in(const Alignment& A);

std::ostream& operator<<(std::ostream& os, const Alignment& a);

