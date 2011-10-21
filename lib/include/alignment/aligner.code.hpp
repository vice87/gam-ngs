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
     
