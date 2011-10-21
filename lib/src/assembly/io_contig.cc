#include "assembly/io_contig.code.hpp"

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
read_contig_in(const std::string& name, const std::string& dir); 

std::istream& 
operator>>(std::istream& is, Contig& ctg);

std::ostream& 
operator<<(std::ostream& os, const Contig& ctg);

const std::vector<Contig>& 
write_fasta(const std::vector<Contig>& ctgs, 
                                const std::string& fasta_filename);

const std::vector<Contig>& 
write_quality(const std::vector<Contig>& ctgs, 
                                const std::string& quality_filename); 

const std::vector<Contig>& 
write_contig(const std::vector<Contig>& ctgs, 
                                const std::string& fasta_filename,
                                const std::string& quality_filename);

const Contig& 
write_fasta(const Contig& ctg, const std::string& fasta_filename);

const Contig& 
write_quality(const Contig& ctg, const std::string& quality_filename); 

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

