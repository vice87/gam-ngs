#ifndef _BANDED_SMITH_WATERMAN_
#define _BANDED_SMITH_WATERMAN_

#include "alignment/my_alignment.hpp"
#include "assembly/contig.hpp"

#define FORCE_MAXGAP_LEN 10
#define DEFAULT_BAND_SIZE 150
#define BSW_MAX_ALIGNMENT 500000

class BandedSmithWaterman
{
    public:
        typedef long int int_type;
        typedef unsigned long int size_type;

    private:
        const ScoreType _match_score;
        const ScoreType _mismatch_score;
        const ScoreType _gap_score;
        const size_type _band_size;

    public:

        BandedSmithWaterman();

        BandedSmithWaterman(
                const ScoreType& match_score,
                const ScoreType& mismatch_score,
                const ScoreType& gap_score,
                const size_type& band_size
                );

        BandedSmithWaterman( const size_type& band_size );

        MyAlignment
        find_alignment(const Contig& a, size_type begin_a, size_type end_a,
                const Contig& b, size_type begin_b, size_type end_b,
				bool force_start = false, bool force_end = false ) const;
};

#endif // _BANDED_SMITH_WATERMAN_

