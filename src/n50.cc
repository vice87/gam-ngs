/* 
 * File:   n50.cc
 * Author: riccardo
 *
 * Created on 1 agosto 2011, 22.39
 */

#include <cstdlib>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>

#define BUFFER_LEN 10240

using namespace std;

long int getNextSequenceLength(istream & fasta)
{
    char b[BUFFER_LEN+1];
    b[BUFFER_LEN] = '\0';
    stringstream buf;
    char c = fasta.peek();
    long int length = 0;
    
    while(!fasta.eof() and (c == ' ' or c == '\n')) 
    {
        fasta.ignore(1);
        if( !fasta.eof() ) c = fasta.peek();
    }
    
    if (!fasta.eof() and c != '>') 
    {
        stringstream ss;
        ss << "next character is " << c;
        return -1;
    }
    
    fasta.getline(b, BUFFER_LEN);
    
    string temp;
    
    while (!fasta.eof() and (fasta.peek() != '>')) 
    {
        fasta >> temp;
        
        for (string::iterator iter = temp.begin(); iter != temp.end(); iter++) 
        {
           if ( (*iter >= 'a' and *iter <= 'z') or (*iter >= 'A' and *iter <= 'Z')) 
               length++;
           else
               return -1;
	}
        
        c = fasta.peek();
        
        while (!fasta.eof() and (c == ' ' or c == '\n')) 
        {
            fasta.ignore(1);
            if (!fasta.eof()) c = fasta.peek();
        }
    }
    
    return length;
}

char * getBaseName( char *path )
{
    char *ptr = strrchr(path, '/');
    return ptr ? ptr + 1 : (char *)path;
}


/*
 * 
 */
int main(int argc, char** argv) 
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <FASTA file>" << std::endl;
        return 1;
    }
    
    ifstream ifs ( argv[1] , ifstream::in );
    long int length;
    long int tot_length;
    
    list< long int > seqLenList;
    
    while( !ifs.eof() )
    {
        length = getNextSequenceLength( ifs );
        seqLenList.push_back( length );
        
        if( length < 0 )
        {
            std::cout << getBaseName(argv[1]) << ": Incorrect FASTA file" << std::endl;
            return 1;
        }
    }
    
    // sort contig lengths
    seqLenList.sort();
    
    // compute total length of the assembly
    tot_length = 0;
    list< long int >::reverse_iterator rit;
    for( rit = seqLenList.rbegin() ; rit != seqLenList.rend(); ++rit ) tot_length += *rit;
    
    // compute N50 and L50
    long int N50 = 0, L50 = 0;
    length = 0;
    
    long int max_length = (seqLenList.rbegin() != seqLenList.rend()) ? *(seqLenList.rbegin()) : 0;
    
    for( rit = seqLenList.rbegin() ; rit != seqLenList.rend(); ++rit )
    {
        if( 2*length >= tot_length ) break;
        length += *rit;
        N50++;
        L50 = *rit;
    }
    
    cout << getBaseName(argv[1]) << " statistics:" << endl 
         << "Total length = " << tot_length 
            << "\tAverage = " << (tot_length * 1.0) / seqLenList.size() 
            << "\tMax = " << max_length << endl
         << "Sequences = " << seqLenList.size() << endl
         << "N50 = " << N50 << endl 
         << "L50 = " << L50 << endl;
    
    ifs.close();
    
    return 0;
}

