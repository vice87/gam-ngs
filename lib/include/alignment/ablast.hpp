#ifndef _ABLAST_
#define _ABLAST_

#define ABLAST_DEFAULT_WORD_SIZE 20

#include <list>
#include <map>
#include <vector>

#include "assembly/contig.hpp"

class ABlast
{

private:
    
    typedef std::map< size_t, std::list<size_t> > HashType;
    typedef std::vector< uint64_t > FoundVectorType;
    
    size_t _word_size; //!< The used word size.
     
    inline 
    size_t sequence_code(const Contig& a, size_t pos)
    {
        size_t code = 0;
        for( size_t i=pos; i < pos + _word_size; i++) code=((LAST_BASE-1)*code)+(a.at(i)).base();
        
        return code;
    }

    inline
    HashType
    build_hash( const Contig& a, uint64_t start, uint64_t end ) 
    {
        HashType a_hash;
        for( size_t i=start; i <= end-_word_size+1; i++ ) a_hash[ sequence_code(a,i) ].push_back(i);

        return a_hash;
    }

    inline
    const FoundVectorType&
    mark_found( FoundVectorType& f_vector, size_t idx_a, size_t idx_b ) 
    {
        if( idx_a < idx_b ) return f_vector;
        
        f_vector.at( idx_a-idx_b ) += 1;
        return f_vector;
    }

    inline
    FoundVectorType
    build_corrispondences_vector(
                    const Contig& a, uint64_t a_start, uint64_t a_end,
                    const Contig& b, uint64_t b_start, uint64_t b_end )
    {
        HashType a_hash = build_hash( a, a_start, a_end );

        FoundVectorType f_vector( a_end-a_start+1, 0 );

        for( size_t b_pos = b_start; b_pos <= b_end-_word_size+1; b_pos++ ) 
        {
            HashType::iterator it = a_hash.find( sequence_code(b,b_pos) );
            
            if( it != a_hash.end() )
            {
                for( std::list<size_t>::const_iterator a_pos = (it->second).begin(); a_pos != (it->second).end(); a_pos++ )
                {
                    size_t idx_a = *a_pos - a_start;
                    size_t idx_b = b_pos - b_start;
                    mark_found( f_vector, idx_a, idx_b);
                }
            }
        }

        return f_vector;
    }

public:

     ABlast();
     
     ABlast( const size_t word_size );

     const ABlast& setWordSize(const size_t word_size);

     size_t getWordSize() const;

     std::list< uint32_t > 
     findHits(const Contig& a, uint64_t a_start, uint64_t a_end, const Contig& b, uint64_t b_start, uint64_t b_end);
};
     
#endif // _ABLAST_
     
