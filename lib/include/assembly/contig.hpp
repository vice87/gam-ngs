#ifndef CONTIGHPP_
#define CONTIGHPP_

#include<string>
#include<vector>

#include "assembly/nucleotide.hpp"

typedef unsigned short QualType;
typedef std::vector<Nucleotide> SeqType;
typedef std::vector<QualType> QualSeqType;

class Contig;

std::istream& operator>>(std::istream&, Contig&);

Contig read_contig(const std::string&, const std::string&, 
                                        const std::string&);

class Contig {
 private:
  std::string _name;
  SeqType _sequence;
  //QualSeqType _quality;

 protected:
  

 public:
  Contig();

  Contig(const Contig& orig);

  Contig(const std::string& name);

  Contig(const std::string& name, const SeqType& sequence,
                                  const QualSeqType& quality);

  Contig(const std::string& name, const SeqType& sequence);

  Contig(const std::string& name, const size_t& size);

  Contig(const size_t& size);

  const std::string &name() const; 

  const Contig& set_name(const std::string& new_name); 

  size_t size() const;

  const Contig& operator=(const Contig& orig); 

  bool operator==(const Contig& ctg) const;

  bool operator!=(const Contig& ctg) const;

  const Nucleotide& operator[](const size_t& index) const;

  Nucleotide& operator[](const size_t& index);

  const Nucleotide& at(const size_t& index) const;

  Nucleotide& at(const size_t& index);

  //const QualType& qual(const size_t& index) const;

  //QualType& qual(const size_t& index);

  SeqType sequence(const size_t& index) const;

  SeqType sequence(const size_t& index, const size_t& length) const;
  
  size_t resize(const size_t& size);

  friend std::istream& operator>>(std::istream&, Contig&);
  friend Contig read_contig(const std::string&, const std::string&, 
                                                        const std::string&);
};

Contig chop_begin(const Contig& ctg, 
                    const size_t& preserve_from);

Contig chop_end(const Contig& ctg, 
                    const size_t& preserve_until);

Contig chop_borders(const Contig& ctg, 
                    const size_t& preserve_from,
                    const size_t& preserve_until);

Contig reverse_complement(const Contig& ctg); 

class Alignment;

Contig
merge(const Contig& a, const Contig& b, 
              const Alignment& alignment);

#endif // CONTIG_
