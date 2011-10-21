#include "assembly/contig.code.hpp"

class Contig;

std::istream& operator>>(std::istream&, Contig&);

Contig read_contig(const std::string&, const std::string&, 
                                        const std::string&);

Contig reverse_complement(const Contig& ctg); 

Contig chop_begin(const Contig& ctg, 
                    const size_t& preserve_from);

Contig chop_end(const Contig& ctg, 
                    const size_t& preserve_until);

Contig chop_borders(const Contig& ctg, 
                    const size_t& preserve_from,
                    const size_t& preserve_until);

Contig
merge(const Contig& a, const Contig& b, 
              const Alignment& alignment);

