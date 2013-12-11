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

#ifndef _IO_CONTIG_
#define _IO_CONTIG_

#include <map>
#include <string>

#include "assembly/contig.hpp"
#include "assembly/RefSequence.hpp"

#define SEQ_LINE_LENGTH 60
#define QUAL_LINE_LENGTH 20

QualSeqType
read_quality(const std::string& filename, const std::string& name);

SeqType
read_fasta(const std::string& filename, const std::string& name);

Contig
read_contig(const std::string& filename,
            const std::string& name);

Contig
read_contig(const std::string& fasta_filename,
            const std::string& quality_filename, const std::string& name);

Contig
read_contig_in(const std::string& ctg_name, const std::string& dir);

std::istream&
operator>>(std::istream& is, Contig& ctg);

std::ostream&
operator<<(std::ostream& os, const Contig& ctg);

const std::vector<Contig>&
write_fasta(const std::vector<Contig>& ctgs,
                                const std::string& fasta_filename);

//const std::vector<Contig>&
//write_quality(const std::vector<Contig>& ctgs,
//                                const std::string& quality_filename);

const std::vector<Contig>&
write_contig(const std::vector<Contig>& ctgs,
                                const std::string& fasta_filename,
                                const std::string& quality_filename);

const Contig&
write_fasta(const Contig& ctg, const std::string& fasta_filename);

//const Contig&
//write_quality(const Contig& ctg, const std::string& quality_filename);

const Contig&
write_contig(const Contig& ctg, const std::string& fasta_filename,
                                const std::string& quality_filename);

const Contig&
write_contig(const Contig& ctg, const std::string& fasta_filename);

const Contig&
write_contig_in(const Contig& ctg, const std::string& dir);

void
delete_file_in(const std::string& filename, const std::string& dir);

void
delete_contig_in(const Contig& ctg, const std::string& dir);

void
split_multifasta_sequence_in(const std::string& dir,
                    const std::string& multifasta_filename);

void
merge_multifasta_sequences_in(const std::string& dir,
                    const std::string& multifasta_filename);

void
readNextContigID( std::istream &is, std::string &ctg_id );

void
readNextSequence( std::istream &is, Contig &ctg );

void
loadSequences( const std::string &file, RefSequence &refSequence, const std::map< std::string, int32_t > &ctg2Id );

#endif // _IO_CONTIG_
