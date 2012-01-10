/* 
 * File:   new-alignments.cc
 * Author: vice
 *
 * Created on 22 dicembre 2011, 15.52
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>

#include "alignment/banded_smith_waterman.hpp"
#include "assembly/io_contig.hpp"
#include "alignment/full_smith_waterman.hpp"
#include "alignment/ablast_new.hpp"

#define BUFFER_LEN 16384

void readNextContigID( std::istream &is, std::string &ctg_id )
{
    char c = is.peek();
    
    while( !is.eof() and (c == ' ' or c == '\n') ) 
    {
        is.ignore(1);
        if( !is.eof() ) c = is.peek();
    }
    
    if( !is.eof() and c != '>' ) 
    {
        std::stringstream ss;
        ss << "Found invalid character: " << c;
	throw std::domain_error(ss.str().c_str());
    }
    
    std::string line, id;
    
    // get name    
    std::getline( is, line );
    id = line.substr(1,line.size()-1);
    
    size_t pos = id.find(' ');
    if( pos != std::string::npos ) id = id.substr(0,pos);
    
    ctg_id = id;
}

void readNextSequence( std::istream &is, Contig &ctg )
{
    if( is.eof() ) return;
    
    ctg.resize(0);
    
    char c('\n');
    UIntType idx = 0;
    
    // while we do not reach a new contig
    while( !is.eof() and c != '>' )
    {
        // read a new char
        is.get(c);
        
        if( c != '\n' and c != '>' and c != ' ' and !is.eof() ) 
        {
            // copy read nucleotide into sequence
            ctg.resize( ctg.size()+1 );
            ctg.at(idx) = c;
            idx++;
        }
    }
    
    if( !is.eof() ) is.unget();
}

void loadMap( std::string file, std::map< std::string, Contig > &pool )
{
    std::ifstream ifs( file.c_str(), std::ifstream::in );
    
    char buffer[BUFFER_LEN];
    ifs.rdbuf()->pubsetbuf( buffer, BUFFER_LEN );
    
    while( !ifs.eof() )
    {
        std::string ctg_name;
        readNextContigID( ifs, ctg_name );
        
        Contig ctg;
        readNextSequence( ifs, ctg );
        
        pool[ ctg_name ] = ctg;
    }
    
    ifs.close();
}

/*
 * 
 */
int main(int argc, char** argv) 
{
    std::string fasta1 = argv[1];
    std::string fasta2 = argv[2];
    
    std::map< std::string, Contig > pool1, pool2;
    
    std::cout << "Loading fasta files... " << std::flush;
    loadMap( fasta1, pool1 );
    loadMap( fasta2, pool2 );
    std::cout << "done.\n" << std::endl;
    
    std::map< std::string, Contig >::const_iterator i,j;
    for( i = pool1.begin(); i != pool1.end(); i++ )
    {
        for( j = pool2.begin(); j != pool2.end(); j++ )
        {
            std::cout << j->first << "(" << (j->second).size() << ")" << " vs. " << i->first << "(" << (i->second).size() << ")" << ":" << std::endl;
            
            FullSmithWaterman sw;
            //MyAlignment swa = sw.find_alignment(i->second,0,(i->second).size()-1,j->second,0,(j->second).size()-1);
            //std::cout << swa << std::endl;
            
            ABlast ablast;
            MyAlignment alignment = ablast.find_alignment( reverse_complement(i->second), j->second );
            std::cout << alignment << std::endl;
                    //sw.find_alignment( i->second, 200, (i->second).size()-1, j->second, 0, (j->second).size()-1 );
            
            std::cout << "(Reversed)"  << j->first << "(" << (j->second).size() << ")" << " vs. " << i->first << "(" << (i->second).size() << ")" << ":" << std::endl;
            Contig revC = reverse_complement(j->second);
            
            alignment = ablast.find_alignment( reverse_complement(i->second), revC );
            std::cout << alignment << std::endl;
        }
    }

    return 0;
}

