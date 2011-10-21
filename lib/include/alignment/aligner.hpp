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
     
