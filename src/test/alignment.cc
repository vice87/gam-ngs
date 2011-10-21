/* 
 * File:   n50.cc
 * Author: riccardo
 *
 * Created on 1 agosto 2011, 22.39
 */

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <algorithm>
#include <stdexcept>

#define BUFFER_LEN 10240

#define DEFAULT_MAX_GAPS 300
#define DEFAULT_MAX_SEARCHED_ALIGNMENT 400000

#ifndef MIN_ALIGNMENT
#define MIN_ALIGNMENT 100
#endif

#ifndef MIN_HOMOLOGY
#define MIN_HOMOLOGY 85
#endif

#ifndef MIN_ALIGNMENT_QUOTIENT
#define MIN_ALIGNMENT_QUOTIENT 0.001
#endif

#include "types.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/ablast.hpp"
#include "alignment/smith_waterman.hpp"
#include "pctg/BestPctgCtgAlignment.hpp"
#include "pool/ContigMemPool.hpp"

BestPctgCtgAlignment findBestAlignment( Contig &ctg1, Contig &ctg2 )
{
    ABlast aligner((UIntType)DEFAULT_MAX_SEARCHED_ALIGNMENT, (UIntType)DEFAULT_MAX_GAPS, (UIntType)DEFAULT_MAX_GAPS);
    Contig workingCopy(ctg2), origCtg(ctg2);
    
    // alignment between ctg1 and ctg2
    Alignment af( aligner.find_alignment(ctg1, 0, workingCopy, 0) );
    
    workingCopy = reverse_complement(origCtg);
    
    // alignment between ctg1 and reverse complement of ctg2
    Alignment ar( aligner.find_alignment(ctg1, 0, workingCopy, 0) );
    
    try{
    std::cout << "FORWARD ALIGNMENT:" << std::endl;
    std::cout << "HOMOLOGY = " << af.homology() << std::endl
                               << "LENGTH = " << af.length() << std::endl
			       << "BEGIN CTG1 = " << af.begin_a() << "; BEGIN CTG2 = " << af.begin_b() << std::endl
			       << "LAST_POS_IN = " << last_pos_in(af).first << "   " << last_pos_in(af).second << std::endl
  			       << "CTG_2_SIZE * MIN_ALIGN_QUOT = " << (((double)ctg2.size())*(MIN_ALIGNMENT_QUOTIENT)) << std::endl
                               << "CTG2 POSITION IN CTG1 = " << af.b_position_in_a() << std::endl
                               << "CTG2 END IN CTG1 = " << af.end_b_in_a() << std::endl
                               << "CTG1 POSITION IN CTG2 = " << af.a_position_in_b() << std::endl
                               << "CTG1 END IN CTG2 = " << af.end_a_in_b() << std::endl
                               << std::endl;
    }catch(...){}

    try{
    std::cout << "REVERSE ALIGNMENT:" << std::endl;
    std::cout << "HOMOLOGY = " << ar.homology() << std::endl
                               << "LENGTH = " << ar.length() << std::endl
			       << "BEGIN CTG1 = " << ar.begin_a() << "; BEGIN CTG2 = " << ar.begin_b() << std::endl
			       << "LAST_POS_IN = " << last_pos_in(ar).first << "   " << last_pos_in(ar).second << std::endl
			       << "CTG_2_SIZE * MIN_ALIGN_QUOT = " << (((double)ctg2.size())*(MIN_ALIGNMENT_QUOTIENT)) << std::endl
                               << "CTG2 POSITION IN CTG1 = " << ar.b_position_in_a() << std::endl
                               << "CTG2 END IN CTG1 = " << ar.end_b_in_a() << std::endl
                               << "CTG1 POSITION IN CTG2 = " << ar.a_position_in_b() << std::endl
                               << "CTG1 END IN CTG2 = " << ar.end_a_in_b() << std::endl
                               << std::endl;
    }catch(...){}

    // return the best alignment found
    if(!ar.is_full())
    {
        if(!af.is_full())
		{
			std::cout << "Both alignments NOT FULL" << std::endl;
			return BestPctgCtgAlignment( Alignment((Contig)ctg1, ctg2, (size_t)0, (size_t)0), false );
		}
		
		std::cout << "FORWARD IS FULL" << std::endl;
        return BestPctgCtgAlignment(af,false);
    }
    else
    {
        if( !af.is_full() || ar.score() > af.score() )
        {
            if( af.is_full() ) std::cout << "FORWARD IS FULL" << std::endl; else std::cout << "FORWARD NOT FULL" << std::endl;
			if( ar.score() > af.score() ) std::cout << "REVERSE HAS BETTER SCORE" << std::endl; else std::cout << "FORWARD HAS BETTER SCORE" << std::endl;
			//ctg = workingCopy;
            return BestPctgCtgAlignment(ar,true);
        }
        
		std::cout << "FORWARD IS FULL AND HAS BETTER SCORE" << std::endl;		
        return BestPctgCtgAlignment(af,false);
    }
}

/*
 * 
 */
int main(int argc, char** argv) 
{
	if(argc != 3) std::cerr << "Usage: " << argv[0] << " <fasta1> <fasta2>" << std::endl;
	
	std::string fasta1 = argv[1];
	std::string fasta2 = argv[2];
	
	ContigMemPool pool1 = ContigMemPool::loadPool(fasta1);
	ContigMemPool pool2 = ContigMemPool::loadPool(fasta2);
	
	ContigMemPool::iterator i,j;
	for( i = pool1.begin(); i != pool1.end(); i++ )
	{
		for( j = pool2.begin(); j != pool2.end(); j++ )
		{
			std::cout << "ALIGNMENT => " << i->first << "(" << (i->second).size() << ") vs. " 
				  << j->first << "(" << (j->second).size() << ")" << std::endl << std::endl;
			
			BestPctgCtgAlignment bestAlign( findBestAlignment( i->second, j->second ) );
			
			if( bestAlign.getAlignment().homology() < MIN_HOMOLOGY ||
				bestAlign.getAlignment().length() < MIN_ALIGNMENT  ||
				bestAlign.getAlignment().length() < ((double)(j->second).size())*MIN_ALIGNMENT_QUOTIENT )
			{
				std::cout << "ALIGNMENT SUCKS!" << std::endl;
			}
			else
			{
				std::cout << "ALIGNMENT IS GOOD!" << std::endl;
				        //<< "HOMOLOGY = " << bestAlign.getAlignment().homology() << std::endl 
					//<< "LENGTH = " << bestAlign.getAlignment().length << std::endl;
					//<< "CTG_2_SIZE * MIN_ALIGN_QUOT" << (((double)(j->second).size())*(MIN_ALIGNMENT_QUOTIENT)) << std::endl;
			}

			std::cout << "HOMOLOGY = " << bestAlign.getAlignment().homology() << std::endl 
				  << "LENGTH = " << bestAlign.getAlignment().length() << std::endl
				  << "CTG2 POSITION IN CTG1 = " << bestAlign.getAlignment().b_position_in_a() << std::endl
				  << "CTG2 END IN CTG1 = " << bestAlign.getAlignment().end_b_in_a() << std::endl
				  << "CTG1 POSITION IN CTG2 = " << bestAlign.getAlignment().a_position_in_b() << std::endl
				  << "CTG1 END IN CTG2 = " << bestAlign.getAlignment().end_a_in_b() << std::endl
				  << std::endl;
		}
	}
}

