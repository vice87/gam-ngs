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

#ifndef CONTIGINPCTGINFO_HPP
#define	CONTIGINPCTGINFO_HPP

#include "types.hpp"
#include "pctg/BestPctgCtgAlignment.hpp"

// interval in ctg structure
typedef struct merge_gap
{
	uint64_t start;
	uint64_t end;
	int64_t gap;
} t_merge_gap;

//! Class storing the properties of a contig inside a paired contig.
class ContigInPctgInfo
{

private:
    int32_t _ctgId;      //!< contig identifier
    UIntType _size;     //!< size of the contig
    int64_t _position;  //!< position inside the paired contig
    bool _reversed;     //!< whether the contig is reverse complemented or not.
    std::list< t_merge_gap > _merge_gaps; //!< extensions/contractions due to merging the contig
    uint64_t _left_cut;  //!< left cut when discarding tail in merging
    uint64_t _right_cut; //!< right cut when discarding tail in merging

public:
    //! A constructor.
    /*!
     * Creates an empty ContigInPctgInfo object.
     */
    ContigInPctgInfo();

    //! A constructor.
    /*!
     * Creates a ContigInPctgInfo object, given a contig ID, its size and a position.
     * \param ctgId     contig identifier
     * \param size      size of the contig
     * \param position  position of the contig inside a paired contig.
     */
    ContigInPctgInfo(const int32_t ctgId, const UIntType &size, const UIntType &position);

    //! A copy constructor.
    /*!
     * Creates a copy of a given ContigInPctgInfo object.
     * \param orig a ContigInPctgInfo object.
     */
    ContigInPctgInfo(const ContigInPctgInfo &orig);

    //! A constructor.
    /*!
     * Creates a ContigInPctgInfo object, given a contig ID and a BestPctgCtgAlignment object.
     * \param ctgId contig's identifier
     * \param bestAlign a BestPctgCtgAlignment object
     *
     * \sa BestPctgCtgAlignment
     */
    ContigInPctgInfo(const int32_t ctgId, const BestPctgCtgAlignment &bestAlign);

    //! Gets sequence identifier.
    /*!
     * \return contig ID.
     */
    int32_t getId() const;

	uint64_t getSize() const;

    //! Gets the position of the first nucleotide of the contig inside the paired contig.
    /*!
     * \return position of the first nucleotide.
     */
    int64_t getFirstNucleotidePos() const;

    //! Gets the position of the last nucleotide of the contig inside the paired contig.
    /*!
     * \return position of the last nucleotide.
     */
    int64_t getLastNucleotidePos() const;

    //! Tells whether the contig is reverse complemented or not.
    /*!
     * \return whether the contig is reverse complemented or not.
     */
    const bool& isReversed() const;

	uint64_t getLeftCut() const;
	uint64_t getRightCut() const;

    //! Sets the position.
    /*!
     * \param pos position of the contig inside a paired contig.
     */
    void setPosition(const UIntType& pos);

	void setLeftCut( const uint64_t &len );
	void setRightCut( const uint64_t &len );

	const std::list< t_merge_gap >& merge_gaps() const;
	std::list< t_merge_gap >& merge_gaps();

	void addMergeGap( uint64_t start, uint64_t end, int64_t gap );
};


#endif	/* CONTIGINPCTGINFO_HPP */

