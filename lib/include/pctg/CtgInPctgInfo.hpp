/*!
 * \file ContigInPctgInfo.hpp
 * \brief Definition of ContigInPctgInfo
 * \details This file contains the definition of the class representing
 * information regarding a contig inside a paired contig.
 */

#ifndef CTGINPCTGINFO_HPP
#define CTGINPCTGINFO_HPP

#include "types.hpp"

//! Class storing the properties of a contig inside a paired contig.
class CtgInPctgInfo
{

private:
    int32_t _ctgId;      //!< contig identifier
    int64_t _start;  //!< position inside the paired contig
    int64_t _end;  //!< position inside the paired contig
    bool _reversed;     //!< whether the contig is reverse complemented or not.
    bool _isMaster;

public:

	CtgInPctgInfo();

	CtgInPctgInfo(const int32_t ctgId, const int64_t start, const int64_t end, const bool reversed, const bool isMaster );

	CtgInPctgInfo(const CtgInPctgInfo &orig);

    int32_t getId() const;
	int64_t getStart() const;
	int64_t getEnd() const;

	bool isReversed() const;
	bool isMaster() const;
};


#endif	/* CTGINPCTGINFO_HPP */

