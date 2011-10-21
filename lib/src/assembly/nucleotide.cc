#include "assembly/nucleotide.code.hpp"

class Nucleotide;

Nucleotide complement(const Nucleotide &orig);

std::ostream& operator<<(std::ostream& os, const Nucleotide& base); 

