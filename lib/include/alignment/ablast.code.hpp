#ifndef _ABLAST_CODE_
#define _ABLAST_CODE_

#include <list>
#include <vector>
#include <map>

#include "alignment/ablast.hpp"
#include "alignment/smith_waterman.hpp"

typedef std::map<ABlast::size_type, 
                 std::list<ABlast::size_type> > HashType;

typedef std::vector< std::pair<ABlast::size_type,
                               ABlast::size_type> > FoundVectorType;

ABlast::ABlast(): Aligner(), _word_size(DEFAULT_WORD_SIZE) {}

ABlast::ABlast(const ABlast::size_type& max_a_gaps,
               const ABlast::size_type& max_b_gaps):
                  Aligner(max_a_gaps,max_b_gaps),
                  _word_size(DEFAULT_WORD_SIZE) {}

ABlast::ABlast(const ScoreType& gap_score,    
               const ABlast::size_type& max_a_gaps,
               const ABlast::size_type& max_b_gaps):
                     Aligner(gap_score,max_a_gaps,max_b_gaps),
                     _word_size(DEFAULT_WORD_SIZE) {}

ABlast::ABlast(const ABlast::size_type& max_alignment): 
                  Aligner(max_alignment),
                  _word_size(DEFAULT_WORD_SIZE) {}

ABlast::ABlast(const ABlast::size_type& max_alignment, 
       const ABlast::size_type& max_a_gaps,
       const ABlast::size_type& max_b_gaps):
                  Aligner(max_alignment,max_a_gaps,max_b_gaps),
                  _word_size(DEFAULT_WORD_SIZE) {}

ABlast::ABlast(const ABlast::size_type& max_alignment, 
       const ScoreType& gap_score,    
       const ABlast::size_type& max_a_gaps,
       const ABlast::size_type& max_b_gaps):
                  Aligner(max_alignment,gap_score,max_a_gaps,max_b_gaps),
                  _word_size(DEFAULT_WORD_SIZE) {}

const ABlast::size_type&
ABlast::word_size() const 
{
  return _word_size;
}

const ABlast&
ABlast::set_word_size(const ABlast::size_type& word_size)
{
  _word_size=word_size;

  return *this;
}

inline
ABlast::size_type
extended_word_size(const ABlast::size_type& word_size)
{
  return 2*word_size;
}

inline
ABlast::size_type
sequence_code(const Contig& a, 
              const ABlast::size_type& begin_pos,
              const ABlast::size_type& word_size)
{
  ABlast::size_type code=0;
  
  for (ABlast::size_type i=begin_pos; i<begin_pos+word_size; i++) {
    code=((LAST_BASE-1)*code)+(a.at(i)).base();
  }

  return code;
}

inline
HashType
build_hash(const Contig& a, const ABlast::size_type& word_size) 
{
  HashType a_hash;

  for (ABlast::size_type i=0; i<a.size()-word_size; i++) {
    a_hash[sequence_code(a,i,word_size)].push_back(i);
  }

  return a_hash;
}

inline
const FoundVectorType&
mark_found(FoundVectorType& f_vector, 
           ABlast::size_type idx_a, const ABlast::size_type& idx_b, 
                       const ABlast::size_type& size_b, 
                       const ABlast::size_type& word_size) 
{

  ABlast::size_type min=idx_a-idx_b+size_b;
  ABlast::size_type max=std::min(min+extended_word_size(word_size),
                                 (ABlast::size_type)f_vector.size());

  for (ABlast::size_type i=min; i<max; i++) {
    f_vector.at(i).first++;
    f_vector.at(i).second=std::min(f_vector.at(i).second,idx_b);
  }

  return f_vector;
}

inline
const FoundVectorType&
mark_found(FoundVectorType& f_vector, 
           const std::list<ABlast::size_type>& idx_a_list,
           const ABlast::size_type& idx_b,
           const ABlast::size_type& size_b, 
           const ABlast::size_type& word_size) 
{
  for (std::list<ABlast::size_type>::const_iterator 
                                 it=idx_a_list.begin();
                                 it!=idx_a_list.end(); 
                                                it++) {
    mark_found(f_vector,*it,idx_b,size_b,word_size);
  }

  return f_vector;
}

inline
FoundVectorType
build_corrispondences_vector(
                  const Contig& a, const Contig& b, 
                   const ABlast::size_type& word_size)
{
  HashType a_hash=build_hash(a,word_size);

  FoundVectorType f_vector(a.size()+b.size(),
           std::pair<ABlast::size_type,ABlast::size_type>(0,b.size()));

  for (ABlast::size_type i=0; i<b.size()-word_size; i++) {
    HashType::iterator it=
            a_hash.find(sequence_code(b,i,word_size));

    if (it!=a_hash.end()) {
      mark_found(f_vector,it->second,i,b.size(),word_size);
    }
  }

  return f_vector;
}

Alignment
ABlast::apply(const Contig& a, const Contig& b, 
                               const bool& b_rev) const
{
  std::vector< std::pair<size_type,size_type> >f_vector=
            build_corrispondences_vector(a, b, _word_size);

  size_type max_begin(0);

  for (size_type i=0; i< f_vector.size(); i++) {
    if (f_vector.at(i).first>f_vector.at(max_begin).first) {
      max_begin=i;
    }
  }

  size_type begin_b=f_vector.at(max_begin).second;
  size_type begin_a=std::max((long unsigned int)0,max_begin+begin_b-b.size());

  SmithWaterman sw(_gap_score,extended_word_size(_word_size),
                           extended_word_size(_word_size)); 

  return sw.apply(a, begin_a+1, b, begin_b+1, b_rev);
}

Alignment
ABlast::apply(const Contig& a, const Contig& b) const
{
  Alignment align(apply(a,b,false)), 
            align_r(apply(a,reverse_complement(b),true));

  if (align.score()>align_r.score()) {
    return align;
  }

  return align_r;
}

Alignment
ABlast::find_alignment(const Contig& a, const size_type& begin_a, 
                    const Contig& b, const size_type& begin_b) const
{
  return apply(a, b,false); 
}

#endif // _ABLAST_CODE_   
