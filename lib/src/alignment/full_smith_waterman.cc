#include <list>
#include <stdexcept>
#include <iostream>

#include "alignment/full_smith_waterman.hpp"

FullSmithWaterman::FullSmithWaterman() : 
        _match_score(MATCH_SCORE), 
        _mismatch_score(MISMATCH_SCORE), 
        _gap_score(GAP_SCORE)
{}

FullSmithWaterman::FullSmithWaterman(
        const ScoreType& match_score, 
        const ScoreType& mismatch_score, 
        const ScoreType& gap_score) :
        _match_score(match_score), _mismatch_score(mismatch_score),
        _gap_score(gap_score)
{}

MyAlignment
FullSmithWaterman::find_alignment(
        const Contig& a, 
        size_type begin_a, 
        size_type end_a, 
        const Contig& b, 
        size_type begin_b, 
        size_type end_b ) const
{
    size_type x_size = end_b - begin_b + 2;
    size_type y_size = end_a - begin_a + 2;
                
    ScoreType* sw_matrix[x_size];
    for( size_type i = 0; i < x_size; i++ ) sw_matrix[i] = new ScoreType[y_size];
    
    for( size_type i = 0; i < x_size; i++ ) sw_matrix[i][0] = 0;
    for( size_type j = 0; j < y_size; j++ ) sw_matrix[0][j] = 0;
    
    // fill SmithWaterman matrix
    for( size_type i = 1; i < x_size; i++ )
    {
        for( size_type j = 1; j < y_size; j++ )
        {
            ScoreType diag_score = sw_matrix[i-1][j-1] + (a.at(begin_a + j-1) == b.at(begin_b + i-1) ? _match_score : _mismatch_score );
            ScoreType up_score = sw_matrix[i-1][j] + _gap_score;
            ScoreType left_score = sw_matrix[i][j-1] + _gap_score;
            
            sw_matrix[i][j] = std::max( std::max(diag_score,up_score), left_score );
        }
    }
    
    // find max score
    size_type max_i = x_size-1;
    size_type max_j = y_size-1;
    ScoreType max_score = sw_matrix[max_i][max_j];
    
    for( size_type j = 1; j < y_size; j++ )
    {
        if( sw_matrix[x_size-1][j] > max_score )
        {
            max_score = sw_matrix[x_size-1][j];
            max_i = x_size-1;
            max_j = j;
        }
    }
        
    for( size_type i = 1; i < x_size; i++ )
    {
        if( sw_matrix[i][y_size-1] > max_score )
        {
            max_score = sw_matrix[i][y_size-1];
            max_i = i;
            max_j = y_size-1;
        }
    }
    
    std::list< AlignmentAlphabet > edit_string;
    
    // traceback to find alignment
    size_type i = max_i;
    size_type j = max_j;
    while( i > 0 && j > 0 )
    {
        ScoreType diag_score = sw_matrix[i-1][j-1] + (a.at(begin_a + j-1) == b.at(begin_b + i-1) ? _match_score : _mismatch_score );
        ScoreType up_score = sw_matrix[i-1][j] + _gap_score;
        //ScoreType left_score = sw_matrix[i][j-1] + _gap_score;
        
        if( sw_matrix[i][j] == diag_score )
        {
            edit_string.push_front( a.at(begin_a + j-1) == b.at(begin_b + i-1) ? MATCH : MISMATCH );
            i--;
            j--;
        }
        else if( sw_matrix[i][j] == up_score )
        {
            edit_string.push_front(GAP_A);
            i--;
        }
        else // sw_matrix[i][j] == left_score
        {
            edit_string.push_front(GAP_B);
            j--;
        }
    }
    
    for( size_type i=0; i < x_size; i++ ) delete[] sw_matrix[i];
    
    MyAlignment sw_alignment( j, i, a.size(), b.size(), max_score, edit_string );
    return sw_alignment;
}