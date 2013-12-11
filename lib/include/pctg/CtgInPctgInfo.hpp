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

