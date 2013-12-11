/*
 *  This file is part of GAM-NGS.
 *  Copyright (c) 2011 by Riccardo Vicedomini <rvicedomini@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Simone Scalabrin <scalabrin@appliedgenomics.org>,
 *  Lars Arverstad <lars.arvestad@scilifelab.se>,
 *  Alberto Policriti <policriti@appliedgenomics.org>,
 *  Alberto Casagrande <casagrande@appliedgenomics.org>
 *
 *  GAM-NGS is an evolution of a previous work (GAM) done by Alberto Casagrande,
 *  Cristian Del Fabbro, Simone Scalabrin, and Alberto Policriti.
 *  In particular, GAM-NGS has been adapted to work on NGS data sets and it has
 *  been written using GAM's software as starting point. Thus, it shares part of
 *  GAM's source code.
 *
 *  GAM-NGS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GAM-NGS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GAM-NGS.  If not, see <http://www.gnu.org/licenses/>.
 *
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

    fasta.getline(b, BUFFER_LEN); // identifier line

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
        std::cerr << "Usage: " << getBaseName(argv[0]) << " <fasta file>" << std::endl;
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
            std::cout << getBaseName(argv[1]) << ": Incorrect fasta file" << std::endl;
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
	long int min_length = (seqLenList.begin() != seqLenList.end()) ? *(seqLenList.begin()) : 0;

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
            << "\tMax = " << max_length
            << "\tMin = " << min_length << endl
         << "Sequences = " << seqLenList.size() << endl
         << "N50 = " << N50 << endl
         << "L50 = " << L50 << endl;

    ifs.close();

    return 0;
}

