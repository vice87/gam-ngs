#ifndef REFSEQUENCE_HPP_
#define REFSEQUENCE_HPP_

#include "assembly/contig.hpp"

struct reference_t {
	std::string RefName;
	int32_t RefLength;
	Contig *Sequence;
};

typedef std::vector< reference_t > RefSequence;

#endif //REFSEQUENCE_HPP_
