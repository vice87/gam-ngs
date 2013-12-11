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

#ifndef _ABLAST_NEW_
#define _ABLAST_NEW_

/*! \file
 *  \brief Definition of ABlast class.
 *  \details This file contains the defintion of a BLAST alike aligner
 *  called ABlast.
 */

#include "my_alignment.hpp"

#define DEFAULT_WORD_SIZE 20

//! A BLAST alike aligner.
/*!
  To align two nucleotide sequences, <i>a</i> and <i>b</i>, ABlast builds the
  hash table of all the words of size <i>word_size</i> of <i>a</i> and
  finds the maximum sequence of that words which are consecutive in <i>b</i>.

  \todo Tranform this class into a template over:
        <ul>
          <li>the used integer class</li>
          <li>a total ordering for alignment</li>
        </ul>

*/
class ABlast
{

public:
    typedef long int int_type;
    typedef unsigned long int size_type;

private:
    const ScoreType _match_score;
    const ScoreType _mismatch_score;
    const ScoreType _gap_score;
    const size_type _band_size;
    size_type _word_size; //!< The used word size.

public:

     //! A costructor.
     /*!
       This costructor initializes the ABlast aligner setting the word size
       to <code>DEFAULT_WORD_SIZE</code>.
      */
     ABlast();

     ABlast(
             const ScoreType& match_score,
             const ScoreType& mismatch_score,
             const ScoreType& gap_score,
             const unsigned int& band_size
     );

     //! A distructor.
     ~ABlast() {}

     //! Set the used word size
     /*!
       This method allows to set the word size used during hash computation.
       @param word_size The new word size.
       @return A references to the current object with the new word size.
      */
     const ABlast& set_word_size(const size_type& word_size);

     //! Get the used word size
     /*!
       This method allows to query the word size used during hash computation.
       @return The word size used during has computation.
      */
     const size_type& word_size() const;

     //! Apply ABlast aligner
     /*!
       Apply the current ABlast aligner to two Contig, <i>a</i> and
       <i>b</i>. If the boolean value of the third parameter is <i>true</i>,
       then <i>b</i> has been complemented and reversed with respect to
       its initial value. Otherwise, <i>b</i> is in its original form.
       @param a A contig to be aligned.
       @param b A contig to be aligned.
       @param b_rev Tell whether <i>b</i> has been reversed or not.
       @return The ``best'' alignment between <i>a</i> and <i>b</i>.
      */
     MyAlignment
     apply(const Contig& a, const Contig& b, const bool& b_rev) const;

     //! Find the ``best'' alignment.
     /*!
       Find the ``best'' alignment between contigs <i>a</i> and <i>b</i>.
       @param a A contig to be aligned.
       @param b A contig to be aligned.
       @return The ``best'' alignment between <i>a</i> and <i>b</i>.
      */
     MyAlignment
     find_alignment( const Contig& a, const Contig& b ) const;
};

#endif // _ABLAST_NEW_

