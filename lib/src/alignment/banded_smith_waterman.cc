#include <list>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include "alignment/banded_smith_waterman.hpp"

pthread_mutex_t mutex_new = PTHREAD_MUTEX_INITIALIZER;

BandedSmithWaterman::BandedSmithWaterman() : 
        _match_score(MATCH_SCORE), 
        _mismatch_score(MISMATCH_SCORE), 
        _gap_score(GAP_SCORE),
        _band_size(DEFAULT_BAND_SIZE)
{}

BandedSmithWaterman::BandedSmithWaterman(
        const ScoreType& match_score, 
        const ScoreType& mismatch_score, 
        const ScoreType& gap_score,
        const size_type& band_size) :
        _match_score(match_score), 
        _mismatch_score(mismatch_score),
        _gap_score(gap_score), 
        _band_size(band_size)
{}

BandedSmithWaterman::BandedSmithWaterman( const size_type& band_size ) :
        _match_score(MATCH_SCORE), 
        _mismatch_score(MISMATCH_SCORE),
        _gap_score(GAP_SCORE), 
        _band_size(band_size)
{}

MyAlignment
BandedSmithWaterman::find_alignment(
        const Contig& a, 
        size_type begin_a, 
        size_type end_a, 
        const Contig& b, 
        size_type begin_b, 
        size_type end_b ) const
{
    int BLOSUM62[5][5] =
    {
       // A   T   C   G   N
       {  8, -4, -4, -4,  0 }, // A
       { -4,  8, -4, -4,  0 }, // T
       { -4, -4,  8, -4,  0 }, // C
       { -4, -4, -4,  8,  0 }, // G
       {  0,  0,  0,  0,  8 }, // N
    };
    
    size_type x_size = end_b - begin_b + 1;
    size_type y_size = (2 * this->_band_size) + 1;
    
    // allocate smith waterman matrix
    ScoreType* sw[x_size];
    
    pthread_mutex_lock(&mutex_new);
    for( size_type i = 0; i < x_size; i++ ) sw[i] = new ScoreType[y_size];
    pthread_mutex_unlock(&mutex_new);
    
    // initialization of the first row
    for( size_type j = 0; j < y_size; j++ )
    {
        int_type pos = begin_a - this->_band_size + j;
        
        if( pos >= 0 && pos < a.size() )
        {
            ScoreType diag = BLOSUM62[a.at(pos).base()][b.at(begin_b).base()]; //(a.at(pos) == b.at(begin_b)) ? this->_match_score : this->_mismatch_score;
            ScoreType up = this->_gap_score;
            ScoreType left = (pos > 0 && j > 0) ? sw[0][j-1] : this->_gap_score;
            
            sw[0][j] = (pos > 0 && j > 0) ? std::max(std::max(diag,up),left) : std::max(up,diag);
        }
    }
    
    // fill SmithWaterman matrix
    for( size_type i = 1; i < x_size; i++ )
    {
        for( size_type j = 0; j < y_size; j++ )
        {
            int_type pos = begin_a + i + j - this->_band_size;
            
            if( pos >= 0 && pos < a.size() )
            {
                if( pos == 0 )
                {
                    ScoreType diag = BLOSUM62[a.at(pos).base()][b.at(begin_b+i).base()]; // ((a.at(pos) == b.at(begin_b+i)) ? this->_match_score : this->_mismatch_score);
                    ScoreType up = (j < y_size-1) ? sw[i-1][j+1] + this->_gap_score : this->_gap_score; 
                    ScoreType left = this->_gap_score;
                                        
                    sw[i][j] = (j < y_size-1) ? std::max(std::max(diag,up),left) : std::max(diag,left);
                }
                else
                {
                    ScoreType diag = sw[i-1][j] + BLOSUM62[a.at(pos).base()][b.at(begin_b+i).base()]; //((a.at(pos) == b.at(begin_b+i)) ? this->_match_score : this->_mismatch_score);
                    ScoreType up = (j < y_size-1) ? sw[i-1][j+1] + this->_gap_score : this->_gap_score;
                    ScoreType left = (j > 0) ? sw[i][j-1] + this->_gap_score : this->_gap_score;
                    
                    if( j < y_size-1 && j > 0 ){ sw[i][j] = std::max(std::max(diag,up),left); }
                    else if( j < y_size-1 ){ sw[i][j] = std::max( diag,up ); } // j == 0
                    else if( j > 0 ){ sw[i][j] = std::max( diag,left ); } // j == y_size-1
                    else { sw[i][j] = diag; } // j == 0 AND j == y_size-1 (only when _band_size == 0)
                }
            }
        }
    }
    
    // find max score
    bool found_max = false;
    int_type max_i = 0, max_j = 0;
    ScoreType max_score = 0;
    
    for( size_type j = 0; j < y_size; j++ )
    {
        int_type pos = begin_a + (x_size-1) + j - this->_band_size;
        
        if( pos >= 0 && pos < a.size() ) /* valid score */
        {
            if( !found_max || sw[x_size-1][j] > max_score )
            {
                found_max = true;
                max_i = x_size-1; max_j = j;
                max_score = sw[x_size-1][j];
            }
        }
    }
    
    //int_type j = 2 * this->_band_size;
    int_type i = ( int_type(a.size()) >= (begin_a+_band_size+1) ) ? int_type(a.size()) - (begin_a+_band_size+1) : 0;
    int_type j = ( int_type(a.size()) >= (begin_a+_band_size+1) ) ? 2 * this->_band_size : 2 * this->_band_size - (begin_a+_band_size+1) + a.size();
    for( ; i < x_size && j >= 0; i++ )
    {
        if( !found_max || sw[i][j] > max_score )
        {
            found_max = true;
            max_i = i; max_j = j;
            max_score = sw[i][j];
        }
        
        j--;
    }
    
    if( !found_max ) return MyAlignment(a,b); // this case shouldn't happen
    
    std::list< AlignmentAlphabet > edit_string;
    
    // traceback to find alignment
    int_type x = max_i; 
    int_type y = max_j;
    int_type pos = begin_a + x + y - this->_band_size;
    
    //std::cout << "max_x=" << x << "\tmax_y=" << y << "pos_a=" << pos << std::endl;
    
    while( x >= 0 && y >= 0 && pos >= 0 )
    {
        if( pos == 0 )
        {
            ScoreType diag = BLOSUM62[a.at(pos).base()][b.at(begin_b+x).base()]; // ((a.at(pos) == b.at(begin_b+x)) ? this->_match_score : this->_mismatch_score);
            ScoreType up = (x > 0 && y < y_size-1) ? sw[x-1][y+1] + this->_gap_score : this->_gap_score; 
            ScoreType left = this->_gap_score;
            
            if( sw[x][y] == diag )
            {
                edit_string.push_front( a.at(pos) == b.at(begin_b + x) ? MATCH : MISMATCH );
                x--;
            }
            else if( y == y_size-1 || sw[x][y] == left )
            {
                edit_string.push_front( GAP_B );
                y--;
            }
            else
            {
                edit_string.push_front( GAP_A );
                x--;
                y++;
            }
        }
        else
        {
            ScoreType diag = (x > 0 ? sw[x-1][y] : 0) + BLOSUM62[a.at(pos).base()][b.at(begin_b+x).base()]; //((a.at(pos) == b.at(begin_b + x)) ? this->_match_score : this->_mismatch_score);
            ScoreType up = (x > 0 && y < y_size-1) ? sw[x-1][y+1] + this->_gap_score : this->_gap_score;
            ScoreType left = (y > 0) ? sw[x][y-1] + this->_gap_score : this->_gap_score;

            if( sw[x][y] == diag )
            {
                edit_string.push_front( a.at(pos) == b.at(begin_b+x) ? MATCH : MISMATCH );
                x--;
            }
            else if( y < y_size-1 && y > 0 && sw[x][y] == up )
            {
                edit_string.push_front( GAP_A );
                x--;
                y++;
            }
            else if( y < y_size-1 && y > 0 ) // left
            {
                edit_string.push_front( GAP_B );
                y--;
            }
            else if( y < y_size-1 ) // y == 0 => up
            {
                edit_string.push_front( GAP_A );
                x--;
                y++;
            }
            else // y == y_size-1 => left
            {
                edit_string.push_front( GAP_B );
                y--;
            }
        }
        
        pos = begin_a + x + y - this->_band_size;
    }
    
    pthread_mutex_lock(&mutex_new);
    for( size_type i=0; i < x_size; i++ ) delete[] sw[i];
    pthread_mutex_unlock(&mutex_new);
    
    MyAlignment sw_alignment( a, pos+1, b, begin_b+x+1, max_score, edit_string );
    return sw_alignment;
}