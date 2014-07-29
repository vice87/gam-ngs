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

#ifndef _IO_CONTIG_CODE_
#define _IO_CONTIG_CODE_

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdexcept>

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include "assembly/io_contig.hpp"

#define BUFFER_LEN 16384

QualSeqType
read_quality(const std::string& filename, const std::string& name)
{
  std::string line, begin_contig=">"+name;
  QualSeqType quality;

  // open file
  std::ifstream qualfile ( filename.c_str() , std::ifstream::in );

  // find contig in multifasta quality
  do {
    getline(qualfile, line);

    if (qualfile.eof()) {
      std::stringstream s;
      s << "Quality sequence of \"" << name << "\" not found.";
      throw std::invalid_argument(s.str());
    }
  } while (line!=begin_contig);

  // get first contig line
  getline(qualfile, line);

  // while we do not reach a new contig
  while ((!qualfile.eof())&&(line.substr(0,1) != ">")) {

    // move read qualities into quality
    std::istringstream is(line);
    QualType qual;

    while (is >> qual) {
      quality.push_back(qual);
    }

    // read a new line
    getline(qualfile, line);
  }

  qualfile.close();
  return quality;
}

SeqType
read_sequence(std::istream& fasta)
{
  char c('\n');
  SeqType sequence;

  if (fasta.eof()) {
    return sequence;
  }

  // while we do not reach a new contig
  while ((!fasta.eof())&&(c != '>')) {

    // read a new char
    fasta.get(c);

    if ((c != '\n')&&(c != '>')&&(c != ' ')&&(!fasta.eof())) {  // I do not understand why I need to request (!fasta.eof())
								// however, if I do not check it, contigs contained in fasta
								// which not end with a newline are enlarged of one nucleotide

      // move read nucleotide into sequence
      sequence.resize(sequence.size()+1);
      sequence[sequence.size()-1]=c;
    }
  }

  if (!fasta.eof()) {
    fasta.unget();
  }

  return sequence;
}

SeqType
read_fasta(const std::string& filename, const std::string& name)
{
  std::string line, begin_contig=">"+name;
  SeqType sequence;

  // open file
  std::ifstream fastafile ( filename.c_str() , std::ifstream::in );

  if (!fastafile.is_open()) {
    std::stringstream ss;

    ss<< "Fasta file "<< filename << " cannot be opened.";

    throw std::domain_error(ss.str().c_str());
  }

  // find contig in multifasta
  do {
    getline(fastafile, line);

    if (fastafile.eof()) {
      std::stringstream s;
      s << "Sequence of \"" << name << "\" not found.";
      throw std::invalid_argument(s.str());
    }

  } while (line!=begin_contig);


  /*
  // get first contig line
  getline(fastafile, line);

  // while we do not reach a new contig
  while ((!fastafile.eof())&&(line.substr(0,1) != ">")) {

    // move read bases into sequence
    unsigned int filled_until=sequence.size();
    sequence.resize(sequence.size()+line.size());
    for (unsigned int i=0; i<line.size(); i++) {
      sequence[i+filled_until]=line[i];
    }

    // read a new line
    getline(fastafile, line);
  }
  */

  sequence=read_sequence(fastafile);

  fastafile.close();
  return sequence;
}

Contig
read_contig(const std::string& filename,
            const std::string& name)
{
  SeqType sequence=read_fasta(filename,name);

  return Contig(name,sequence);
}

Contig
read_contig(const std::string& fasta_filename,
            const std::string& quality_filename, const std::string& name)
{
  SeqType sequence=read_fasta(fasta_filename,name);
  QualSeqType quality=read_quality(quality_filename,name);

  return Contig(name,sequence,quality);
}

Contig
read_contig_in(const std::string& name, const std::string& dir)
{
  std::string fasta_name(dir+std::string("/")
                            +name+std::string(".fasta"));
  std::string qual_name(dir+std::string("/")
                           +name+std::string(".qual"));

  //return read_contig(fasta_name, qual_name, name);
  return read_contig(fasta_name, name);
}

std::istream& operator>>(std::istream& is, Contig& ctg)
{
  std::string line;
  ctg._sequence.resize(0);

  // get name
  getline(is, line);
  ctg._name=line.substr(1,line.size()-1);

  size_t pos = ctg._name.find(' ');
  if( pos > 0 ) ctg._name = ctg._name.substr(0,pos);

  /*
  getline(is, line);

  // while we do not reach a new contig
  while ((!is.eof())&&(line.substr(0,1) != ">")) {

    // move read bases into sequence
    unsigned int filled_until=ctg.size();
    ctg._sequence.resize(ctg.size()+line.size());
    for (unsigned int i=0; i<line.size(); i++) {
      ctg._sequence[i+filled_until]=line[i];
    }

    // read a new line
    getline(is, line);
  }
  */

  ctg._sequence=read_sequence(is);
  //ctg._quality.resize( ctg._sequence.size() );

  return is;
}

std::ostream& operator<<(std::ostream& os, const Contig& ctg)
{
  os << ">" << ctg.name();

  size_t i=0;
  while (i<ctg.size()) {
    os << std::endl;
    size_t j=0;
    while ((i<ctg.size())&&(j<SEQ_LINE_LENGTH)) {
      os << ctg.at(i);
      i++;
      j++;
    }
  }

  return os;
}

const std::vector<Contig>&
write_fasta(const std::vector<Contig>& ctgs,
                                const std::string& fasta_filename)
{

  // write sequences
  std::ofstream seqfile ( fasta_filename.c_str() , std::ofstream::out );
  for (unsigned int i=0; i< ctgs.size(); i++) {
    seqfile << ctgs[i] << std::endl;
  }
  seqfile.close();

  return ctgs;
}

//const std::vector<Contig>&
//write_quality(const std::vector<Contig>& ctgs,
//                                const std::string& quality_filename)
//{
//  // write quality
//  std::ofstream qualfile ( quality_filename.c_str() , std::ofstream::out );
//  for (unsigned int i=0; i< ctgs.size(); i++) {
//    qualfile << ">" << ctgs[i].name();
//
//    unsigned int k=0;
//    while (k<ctgs[i].size()) {
//      qualfile << std::endl;
//      unsigned int j=0;
//      while ((k<ctgs[i].size())&&(j<QUAL_LINE_LENGTH)) {
//        qualfile << ctgs[i].qual(k++) << " ";
//        j++;
//      }
//    }
//    qualfile << std::endl;
//  }
//
//  qualfile.close();
//
//  return ctgs;
//}

const std::vector<Contig>&
write_contig(const std::vector<Contig>& ctgs,
                                const std::string& fasta_filename,
                                const std::string& quality_filename)
{
   write_fasta(ctgs,fasta_filename);
   //write_quality(ctgs,quality_filename);

   return ctgs;
}

const Contig&
write_fasta(const Contig& ctg, const std::string& fasta_filename)
{
  std::vector<Contig> v(1);

  v[0]=ctg;

  write_fasta(v, fasta_filename);

  return ctg;
}

const Contig&
write_quality(const Contig& ctg, const std::string& quality_filename)
{
  std::vector<Contig> v(1);

  v[0]=ctg;

  //write_quality(v, quality_filename);

  return ctg;
}

const Contig&
write_contig(const Contig& ctg, const std::string& fasta_filename,
                                const std::string& quality_filename)
{
  std::vector<Contig> v(1);

  v[0]=ctg;

  write_fasta(v, fasta_filename);
  //write_quality(v, quality_filename);

  return ctg;
}

const Contig&
write_contig(const Contig& ctg, const std::string& fasta_filename)
{
  std::vector<Contig> v(1);

  v[0]=ctg;

  write_fasta(v, fasta_filename);

  return ctg;
}

const Contig&
write_contig_in(const Contig& ctg, const std::string& dir)
{
  std::string fasta_name(dir+std::string("/")
                            +ctg.name()+std::string(".fasta"));
  std::string qual_name(dir+std::string("/")
                           +ctg.name()+std::string(".qual"));

  //write_contig(fasta_name, qual_name, ctg.name());
  write_contig(ctg, fasta_name);

  return ctg;
}

void
split_multifasta_sequence_in(const std::string& dir,
                    const std::string& multifasta_filename)
{

#ifdef _VERBOSE_OUTPUT_
  std::cout << "Splitting multifasta file \"" << multifasta_filename
            << "\" in directory \""<< dir << "\"..." <<std::flush;
#endif

  std::string line;

  if ((mkdir(dir.c_str(),S_IRWXU)==-1)&&((errno!=EEXIST)||(errno==ENOTDIR))) {
    throw std::runtime_error("I cannot create the directory");
  }

  // open file
  std::ifstream fastafile ( multifasta_filename.c_str() , std::ifstream::in );

  // get first contig line
  getline(fastafile, line);

  // while we do not reach a new contig
  while (!fastafile.eof()) {

    if (line.substr(0,1) != std::string(">")) {
      throw std::invalid_argument("Wrong multifasta format");
    }

    size_t space_idx=line.find_first_of(" ");

    if (std::string::npos==space_idx) {
      space_idx=line.length();
    }

    std::string ctg_name=line.substr(1,space_idx-1);
    SeqType sequence;

    getline(fastafile, line);
    // while we do not reach a new contig
    while ((!fastafile.eof())&&(line.substr(0,1) != ">")) {

      // move read bases into sequence
      size_t filled_until=sequence.size();
      sequence.resize(sequence.size()+line.size());
      for (unsigned int i=0; i<line.size(); i++) {
        sequence[i+filled_until]=line[i];
      }

      // read a new line
      getline(fastafile, line);
    }

    std::string fastaname(dir+std::string("/")+
                          ctg_name+std::string(".fasta"));

    write_contig(Contig(ctg_name,sequence), fastaname);
  }

  fastafile.close();

#ifdef _VERBOSE_OUTPUT_
  std::cout << "done" << std::endl;
#endif

}

void
merge_multifasta_sequences_in(const std::string& dir,
                    const std::string& multifasta_filename)
{
#ifdef _VERBOSE_OUTPUT_
  std::cout << "Merging sequences contained by directory \"" << dir
            << "\" into multifasta file \"" << multifasta_filename
	    << "\"..." << std::flush;
#endif

  std::ofstream file( multifasta_filename.c_str() );

  DIR *dp;
  struct dirent *dirp;

  if ((dp  = opendir(dir.c_str())) == NULL) {
    throw std::runtime_error("I cannot open the directory");
  }

  while ((dirp = readdir(dp)) != NULL) {
    std::string name(dir+std::string("/")+std::string(dirp->d_name));
    size_t last_dot=name.find_last_of('.');

    if ((last_dot!=std::string::npos)&&
        (name.substr(last_dot)==std::string(".fasta"))) {

      std::ifstream contig_file(name.c_str());
      Contig ctg;

      contig_file >> ctg;

      contig_file.close();

      file << ctg << std::endl;
    }
  }

  closedir(dp);
  file.close();
#ifdef _VERBOSE_OUTPUT_
  std::cout << "done" << std::endl;
#endif

}

void
delete_file_in(const std::string& filename, const std::string& dir)
{
  std::string name(dir+std::string("/")+filename);

  remove(name.c_str());
}

void
delete_contig_in(const Contig& ctg, const std::string& dir)
{
  std::string fasta_name(ctg.name()+std::string(".fasta"));
  std::string qual_name(ctg.name()+std::string(".qual"));

  delete_file_in(fasta_name,dir);
  //delete_file_in(qual_name,dir);
}


void readNextContigID( std::istream &is, std::string &ctg_id )
{
	char c = is.peek();

	while( is.good() and (c == ' ' or c == '\n') )
	{
		is.ignore(1);
		if( is.good() ) c = is.peek();
	}

	if( is.good() and c != '>' )
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

	char c('\n');
	size_t idx = 0;

	// while we do not reach a new contig
	while( !is.eof() and c != '>' )
	{
		// read a new char
		is.get(c);

		if( c != '\n' and c != '>' and c != ' ' and !is.eof() )
		{
			// copy read nucleotide into sequence
			if(idx >= ctg.size()) ctg.resize(idx+1);

		  ctg.at(idx) = c;
		  idx++;
		}
	}

	if( !is.eof() ) is.unget();
}



size_t
loadSequences( const std::string &file, RefSequence &refSequence, const std::map< std::string, int32_t > &ctg2Id )
{
	std::ifstream ifs( file.c_str(), std::ifstream::in );

	char buffer[BUFFER_LEN];
	ifs.rdbuf()->pubsetbuf( buffer, BUFFER_LEN );
	
	size_t num = 0;

	while( !ifs.eof() )
	{
		std::string ctg_name;
		readNextContigID( ifs, ctg_name );

		std::map< std::string, int32_t >::const_iterator it = ctg2Id.find( ctg_name );

		Contig *ctg = new Contig( ctg_name, refSequence[it->second].RefLength );
		readNextSequence( ifs, *ctg );

		refSequence[it->second].Sequence = ctg;
		
		++num;
	}

	ifs.close();
	
	return num;
}

/*
void
delete_contig_in(const Contig& ctg, const std::string& dir)
{
  std::string fasta_name(dir+std::string("/")
                            +ctg.name()+std::string(".fasta"));
  std::string qual_name(dir+std::string("/")
                           +ctg.name()+std::string(".qual"));

  delete_file_in(fasta_name,dir);
  //delete_file_in(qual_name,dir);
}
*/

#endif // _IO_CONTIG_CODE_
