#ifndef _SMITH_WATERMAN_CODE_
#define _SMITH_WATERMAN_CODE_

#include <list>
#include <stdexcept>

#include <stdlib.h>

#include "alignment/smith_waterman.hpp"

void
SmithWaterman::compute_alignment_matrix(
                         const Contig& a, SmithWaterman::size_type& begin_a, 
                         const Contig& b, SmithWaterman::size_type& begin_b, 
                         ScoreType* a_matrix[],
                         const SmithWaterman::size_type& y_size,
                         const SmithWaterman::size_type& x_size) const 
{
    
  for (SmithWaterman::size_type x=0; x< x_size; x++) {
    if (x>=this->_max_b_gaps) {
      SmithWaterman::size_type a_pos=get_a_pos(0,x,begin_a);
      SmithWaterman::size_type b_pos=get_b_pos(0,begin_b);
  
      if ((b_pos<b.size())&&(a_pos<a.size())) {
        a_matrix[0][x]=match_score(a, a_pos, b, b_pos);
      }
    }
  }

  ScoreType m_score;
  for (SmithWaterman::size_type y=1; y < y_size; y++) {
    SmithWaterman::size_type min_vect=0;
   
    if (y<this->_max_b_gaps) {
      min_vect=this->_max_b_gaps-y;
    }

    if (min_vect+1<x_size) {
      SmithWaterman::size_type a_pos=get_a_pos(y,min_vect,begin_a);
      SmithWaterman::size_type b_pos=get_b_pos(y,begin_b);
  
      if ((b_pos<b.size())&&(a_pos<a.size())) {
    
        m_score= match_score(a, a_pos, b, b_pos);

        if ( get_a_pos(y-1,min_vect+1,begin_a)<a.size()) {
          if (min_vect==0) {
            a_matrix[y][min_vect]=std::max(m_score,
                                         a_matrix[y-1][min_vect+1]+
                                                          gap_score());
          } else {
            a_matrix[y][min_vect]=std::max(a_matrix[y-1][min_vect]+m_score,
                                         a_matrix[y-1][min_vect+1]+
                                                         gap_score());
          }
        } else {
          a_matrix[y][min_vect]=m_score;
        }
    
        for (SmithWaterman::size_type x=min_vect+1; x<x_size-1; x++) {
 
          SmithWaterman::size_type a_pos=get_a_pos(y,x,begin_a);
  
          if (a_pos<a.size()) {
            m_score= match_score(a, a_pos, b, b_pos);

            a_matrix[y][x]=std::max(std::max(a_matrix[y-1][x]+m_score,
                                          a_matrix[y-1][x+1]+gap_score()),
                                          a_matrix[y][x-1]+gap_score());
          }
        }

        SmithWaterman::size_type a_pos=get_a_pos(y,x_size-1,begin_a);

        if (a_pos<a.size()) {
          m_score= match_score(a, a_pos, b, b_pos);

          a_matrix[y][x_size-1]=
                         std::max(a_matrix[y-1][x_size-1]+m_score,
                                a_matrix[y][x_size-2]+gap_score());               
        }
      }
    }
  }
}

inline
Alignment 
SmithWaterman::compute_alignment_from_matrix(
                         const Contig& a, SmithWaterman::size_type& begin_a, 
                         const Contig& b, SmithWaterman::size_type& begin_b, 
                         ScoreType* a_matrix[],
                         const SmithWaterman::size_type& y_size,
                         const SmithWaterman::size_type& x_size,
                         const bool& b_rev) const 
{

  SmithWaterman::size_type x=0;
  SmithWaterman::size_type y=0;
  ScoreType max=0;

  for (SmithWaterman::size_type y_s=0; y_s<y_size ; y_s++) {
    for (SmithWaterman::size_type x_s=0; x_s<x_size; x_s++) {
      if ((get_a_pos(y_s,x_s,begin_a)<a.size())&&
          (get_b_pos(y_s,begin_b)<b.size())&&
          (max<a_matrix[y_s][x_s])) {
        max=a_matrix[y_s][x_s];
        x=x_s;
        y=y_s;
      }
    }
  }

  std::list<AlignmentAlphabet> sequence_list;
  std::list<ScoreType> m_score_list;

  bool done=false;
  while ((!((y==0)||(get_a_pos(y,x,begin_a)==begin_a)))&&(!done)) {
    ScoreType m_score=a_matrix[y][x]-a_matrix[y-1][x];

    SmithWaterman::size_type a_pos=get_a_pos(y,x,begin_a);
    SmithWaterman::size_type b_pos=get_b_pos(y,begin_b);

    if (m_score==match_score(a, a_pos, b, b_pos)) {
      if (m_score>0) {
        sequence_list.push_front(MATCH);
      } else {
        sequence_list.push_front(MISMATCH);
      }
      m_score_list.push_front(m_score);
      y--;
    } else { 
     if ((x+1<x_size)&&
        (a_matrix[y-1][x+1]+gap_score()==a_matrix[y][x])) {
       sequence_list.push_front(GAP_A);
       y--;
       x++;
     } else {
       if (a_matrix[y][x-1]+gap_score()==a_matrix[y][x]) {
         sequence_list.push_front(GAP_B);
         x--;
       } else {
          done=true;
        }
      }
    }
  }

  if (!done) {
    if (a_matrix[y][x]>=0) {
      sequence_list.push_front(MATCH);
    } else {
      sequence_list.push_front(MISMATCH);
    }
    m_score_list.push_front(a_matrix[y][x]);
  }

  SmithWaterman::size_type a_pos=get_a_pos(y,x,begin_a);
  SmithWaterman::size_type b_pos=get_b_pos(y,begin_b);

  Alignment max_A(a,b,a_pos,b_pos,b_rev);

  for (std::list<AlignmentAlphabet>::iterator j=sequence_list.begin();
                                           j!=sequence_list.end(); j++) {
    switch (*j) {
     case MATCH:
       max_A.match(m_score_list.front());
       m_score_list.pop_front();
       break;
     case MISMATCH:
       max_A.mismatch(m_score_list.front());
       m_score_list.pop_front();
       break;
     case GAP_A:
       max_A.a_gaps(gap_score());
       break;
     case GAP_B:
       max_A.b_gaps(gap_score());
       break;
    }
  }

  return max_A;
}

SmithWaterman::SmithWaterman() {}
SmithWaterman::SmithWaterman(const SmithWaterman::size_type& max_a_gaps, 
           const SmithWaterman::size_type& max_b_gaps):  
                               Aligner(max_a_gaps,max_b_gaps) {}

SmithWaterman::SmithWaterman(const ScoreType& gap_score, 
                    const SmithWaterman::size_type& max_a_gaps, 
                    const SmithWaterman::size_type& max_b_gaps): 
                          Aligner(gap_score, max_a_gaps, max_b_gaps) {}

SmithWaterman::SmithWaterman(const SmithWaterman::size_type& max_alignment): 
                                              Aligner(max_alignment) {}

SmithWaterman::SmithWaterman(const SmithWaterman::size_type& max_alignment, 
                             const SmithWaterman::size_type& max_a_gaps, 
                             const SmithWaterman::size_type& max_b_gaps):  
                             Aligner(max_alignment,max_a_gaps,max_b_gaps) {}

SmithWaterman::SmithWaterman(const SmithWaterman::size_type& max_alignment, 
                             const ScoreType& gap_score, 
                             const SmithWaterman::size_type& max_a_gaps,  
                             const SmithWaterman::size_type& max_b_gaps): 
                  Aligner(max_alignment, gap_score,max_a_gaps,max_b_gaps) {} 

Alignment
SmithWaterman::apply(const Contig& a, SmithWaterman::size_type begin_a,
               const Contig& b, SmithWaterman::size_type begin_b,
               const bool& b_rev) const
{ 
  if ((a.size()<begin_a)||(b.size()<begin_b)) {
    return Alignment(a,b,begin_a,begin_b,b_rev);
  }

  if (begin_a<this->_max_a_gaps) {
    begin_a=0;
  } else {
    begin_a-=this->_max_a_gaps;
  }

  if (begin_b<this->_max_b_gaps) {
    begin_b=0;
  } else {
    begin_b-=this->_max_b_gaps;
  }

  SmithWaterman::size_type x_size(a.size()-begin_a);
  SmithWaterman::size_type y_size(b.size()-begin_b);

  if (((begin_a==0)&&(begin_b!=0))||
      ((begin_a!=0)&&(begin_b==0))) {
    if (x_size<y_size) {
      y_size=x_size;
    } else {
      x_size=y_size;
    }
  }

  if (this->_max_a_gaps+this->_max_b_gaps!=0) {
    x_size=this->_max_a_gaps+this->_max_b_gaps;
  }

  if ((max_alignment()!=0)&&(y_size>max_alignment())) {
    y_size=max_alignment();
  }

  ScoreType* a_matrix[y_size];

  for (SmithWaterman::size_type y=0; y<y_size; y++) {
    a_matrix[y]=(ScoreType *)malloc(sizeof(ScoreType)*x_size);

    if (a_matrix[y]==NULL) {
      throw std::runtime_error("SmithWaterman function cannot allocate enough memory.");
    }

    for (SmithWaterman::size_type x=0; x< x_size; x++) {
      a_matrix[y][x]=0;
    }
  }

  compute_alignment_matrix(a, begin_a, b, begin_b, 
                               a_matrix, y_size, x_size); 

  Alignment result(compute_alignment_from_matrix(a, 
                               begin_a, b, begin_b,
                                a_matrix, y_size, x_size, b_rev)); 

  for (SmithWaterman::size_type y=0; y<y_size; y++) {
    free(a_matrix[y]);
  }

  return result;
}

Alignment
SmithWaterman::apply(const Contig& a, SmithWaterman::size_type begin_a,
               const Contig& b, SmithWaterman::size_type begin_b) const
{ 
  return apply(a, begin_a, b, begin_b, false);
}

#endif // _SMITH_WATERMAN_CODE_

