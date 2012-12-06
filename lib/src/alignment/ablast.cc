#include <stdint.h>

#include "alignment/ablast.hpp"

ABlast::ABlast() : _word_size(ABLAST_DEFAULT_WORD_SIZE) {}


ABlast::ABlast(const size_t word_size) : _word_size(word_size) {}


std::list< uint32_t >
ABlast::findHits(const Contig& a, uint64_t a_start, uint64_t a_end, const Contig& b, uint64_t b_start, uint64_t b_end)
{
    std::list< uint32_t > hitsList;
    uint64_t max_score(0);
    
    if( a.size() == 0 || b.size() == 0 ) return hitsList;
    
    if( a_end >= a.size() ) a_end = a.size()-1;
    if( b_end >= b.size() ) b_end = b.size()-1;
    
    if( a_start > a_end || b_start > b_end ) return hitsList;
    if( a_end + 1 < _word_size + a_start || b_end + 1 < _word_size + b_start ) return hitsList;
    
    std::vector< uint64_t > f_vector = this->build_corrispondences_vector( a, a_start, a_end, b, b_start, b_end );
    
    // find best hits and fill the output list
    for( size_t i=0; i < f_vector.size(); i++ )
    {
        if( f_vector[i] == 0 ) continue;
        
        if( f_vector[i] > max_score )
        {
            max_score = f_vector[i];
            
            hitsList.clear();
            hitsList.push_back(a_start + i);
        }
        else if( f_vector[i] == max_score )
        {
            hitsList.push_back(a_start + i);
        }
    }
    
    return hitsList;
}

