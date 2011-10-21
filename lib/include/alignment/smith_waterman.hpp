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
     
