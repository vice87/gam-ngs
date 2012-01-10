#ifndef _MY_ALIGNMENT_
#define _MY_ALIGNMENT_

/*! \file
 *  \brief Definition of Alignment class.
 *  \details This file contains the defintion of Alignment.
 */

#include <list>
#include <vector>
#include <iostream>

#include "types.hpp"
#include "assembly/contig.hpp"

#define GAP_SCORE -16
#define MATCH_SCORE 2
#define MISMATCH_SCORE -1

#ifndef FULL_ALIGNMENT_SIZE_PERCENTAGE
#define FULL_ALIGNMENT_SIZE_PERCENTAGE 95.00
#endif

typedef long int ScoreType;

typedef enum {
  GAP_A,
  GAP_B,
  MATCH,
  MISMATCH
} __attribute__((packed)) AlignmentAlphabet;


class MyAlignment 
{
    friend std::ostream& operator<<(std::ostream& os, const MyAlignment& a);
    
public:
    typedef long int int_type;
    typedef unsigned long int size_type;
    typedef std::vector<AlignmentAlphabet> SeqType;

private:
    Contig _a;
    Contig _b;
    size_type _begin_a;
    size_type _begin_b;
    //bool _b_rev;
    SeqType _sequence;
    ScoreType _score;
    
public:
    
    MyAlignment();
    MyAlignment( const MyAlignment& orig );
    MyAlignment( const Contig &a, const Contig &b );
    MyAlignment( const Contig &a, size_type begin_a, const Contig &b, size_type begin_b );
    MyAlignment( const Contig &a, size_type begin_a, const Contig &b, size_type begin_b, ScoreType score, const std::list<AlignmentAlphabet> edit_string );
    
    const Contig& a() const;
    const Contig& b() const;
    
    size_type begin_a() const;
    size_type begin_b() const;
    
    const SeqType& sequence() const;
    
    size_type length() const;
    
    ScoreType score() const;
    
    RealType homology() const; 
    
    int_type b_position_in_a() const;
    int_type a_position_in_b() const;
    
    int_type end_a_in_b() const; 
    int_type end_b_in_a() const;
    
    const MyAlignment& operator=(const MyAlignment& orig);
};

bool first_match_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos );
bool last_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos );
bool last_match_pos( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> &pos ); 
bool gaps_before_last_match( const MyAlignment& A, std::pair<MyAlignment::size_type,MyAlignment::size_type> gaps );

#endif // _MY_ALIGNMENT_
