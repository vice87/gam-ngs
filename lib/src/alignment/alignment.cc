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

