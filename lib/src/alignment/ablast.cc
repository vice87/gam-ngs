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

#include <stdint.h>

#include "alignment/ablast.hpp"

ABlast::ABlast() : _word_size(ABLAST_DEFAULT_WORD_SIZE) {}


ABlast::ABlast(const size_t word_size) : _word_size(word_size) {}


std::list< uint32_t >
ABlast::findHits(const Contig& a, uint64_t a_start, uint64_t a_end, const Contig& b, uint64_t b_start, uint64_t b_end)
{
    std::list< uint32_t > hitsList;
    uint64_t max_score(0);

    if( a.size() == 0 || b.size() == 0 ) return hitsList;

    if( a_end >= a.size() ) a_end = a.size()-1;
    if( b_end >= b.size() ) b_end = b.size()-1;

    if( a_start > a_end || b_start > b_end ) return hitsList;
    if( a_end + 1 < _word_size + a_start || b_end + 1 < _word_size + b_start ) return hitsList;

    std::vector< uint64_t > f_vector = this->build_corrispondences_vector( a, a_start, a_end, b, b_start, b_end );

    // find best hits and fill the output list
    for( size_t i=0; i < f_vector.size(); i++ )
    {
        if( f_vector[i] == 0 ) continue;

        if( f_vector[i] > max_score )
        {
            max_score = f_vector[i];

            hitsList.clear();
            hitsList.push_back(a_start + i);
        }
        else if( f_vector[i] == max_score )
        {
            hitsList.push_back(a_start + i);
        }
    }

    return hitsList;
}

