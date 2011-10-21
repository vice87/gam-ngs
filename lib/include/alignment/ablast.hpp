#ifndef _ABLAST_
#define _ABLAST_

/*! \file
 *  \brief Definition of ABlast class.
 *  \details This file contains the defintion of a BLAST alike aligner
 *  called ABlast.
 */

#define DEFAULT_WORD_SIZE 20

#include "alignment/aligner.hpp"

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
class ABlast : public Aligner
{
  public:
/*     typedef Aligner::int_type int_type;
     typedef Aligner::size_type size_type;
*/
  private:
     size_type _word_size; //!< The used word size.

  public:

     //! A costructor.
     /*!
       This costructor initializes the ABlast aligner setting the word size
       to <code>DEFAULT_WORD_SIZE</code>.
      */
     ABlast();

     //! A costructor.
     /*!
       This costructor has been implemented to get the same API of 
       SmithWaterman. It actually discards the parameters and
       calls ABlast().
       \todo Remove such costructor.
      */
     ABlast(const size_type& max_a_gaps, const size_type& max_b_gaps);  
     
     //! A costructor.
     /*!
       This costructor has been implemented to get the same API of 
       SmithWaterman. It actually discards the parameters and
       calls ABlast().
       \todo Remove such costructor.
      */
     ABlast(const size_type& max_alignment);
     
     //! A costructor.
     /*!
       This costructor has been implemented to get the same API of 
       SmithWaterman. It actually discards the parameters and
       calls ABlast().
       \todo Remove such costructor.
      */
     ABlast(const ScoreType& gap_score, 
                const size_type& max_a_gaps, const size_type& max_b_gaps); 

     //! A costructor.
     /*!
       This costructor has been implemented to get the same API of 
       SmithWaterman. It actually discards the parameters and
       calls ABlast().
       \todo Remove such costructor.
      */
     ABlast(const size_type& max_alignment, const size_type& max_a_gaps, 
                                      const size_type& max_b_gaps);  

     //! A costructor.
     /*!
       This costructor has been implemented to get the same API of 
       SmithWaterman. It actually discards the parameters and
       calls ABlast().
       \todo Remove such costructor.
      */
     ABlast(const size_type& max_alignment, const ScoreType& gap_score, 
                const size_type& max_a_gaps, const size_type& max_b_gaps); 

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
     Alignment
     apply(const Contig& a, const Contig& b, const bool& b_rev) const;

     //! Apply ABlast aligner
     /*! 
       Apply the current ABlast aligner to two Contig, <i>a</i> and 
       <i>b</i>. It return the ``best'' alignment between <i>a</i>, <i>b</i>,
       and the reverse complement of <i>b</i>.
       @param a A contig to be aligned.
       @param b A contig to be aligned.
       @return The ``best'' alignment between <i>a</i>, <i>b</i>, and 
               the reverse complement of <i>b</i>.
      */
     Alignment
     apply(const Contig& a, const Contig& b) const;

     //! Find the ``best'' alignment.
     /*! 
       Find the ``best'' alignment between contigs <i>a</i> and <i>b</i>.
       @param a A contig to be aligned.
       @param b A contig to be aligned.
       @return The ``best'' alignment between <i>a</i> and <i>b</i>. 
      */
     Alignment
     find_alignment(const Contig& a, const size_type& begin_a, 
                    const Contig& b, const size_type& begin_b) const; 
};
     
#endif // _ABLAST_
     
