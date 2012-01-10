#ifndef _FULL_SMITH_WATERMAN_
#define _FULL_SMITH_WATERMAN_

#include "alignment/my_alignment.hpp"
#include "assembly/contig.hpp"

class FullSmithWaterman 
{
    public:
        //typedef Aligner::int_type int_type;
        //typedef Aligner::size_type size_type;
        typedef long int int_type;
        typedef unsigned long int size_type;
        
    private:
        const ScoreType _match_score;
        const ScoreType _mismatch_score;
        const ScoreType _gap_score;
        
    public:
        
        FullSmithWaterman();
        
        FullSmithWaterman(
                const ScoreType& match_score,
                const ScoreType& mismatch_score,
                const ScoreType& gap_score
                );
        
        MyAlignment
        find_alignment(const Contig& a, size_type begin_a, size_type end_a,
                const Contig& b, size_type begin_b, size_type end_b) const;
};
     
#endif // _FULL_SMITH_WATERMAN_
     
