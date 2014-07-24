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
 * \file Read.hpp
 * \brief Definition of Read class.
 * \details This file contains the definition of the class representing a read.
 */

#ifndef READS_HPP
#define	READS_HPP

#include <map>
#include <string>

#include "bam/MultiBamReader.hpp"

#include "types.hpp"
#include "google/sparse_hash_map"

using namespace BamTools;
using google::sparse_hash_map;

//! Class implementing a read.
class Read
{

private:
    int32_t _contigId;   //!< contig's identifier.
    int32_t _startPos;   //!< starting position (0-based) of the read inside the contig.
    int32_t _endPos;     //!< ending position (0-based) of the read inside the contig (half-open interval).
    bool _isRev;        //!< whether the read is reverse complemented or not.

public:
    //! A constructor with no arguments.
    Read();

    //! A copy constructor.
    /*!
     * Creates a copy of a read.
     * \param orig a Read object.
     */
    Read(const Read &orig);

    //! A constructor which sets the read's attributes.
    /*!
     * \param ctg       contig's identifier
     * \param sPos      starting position (0-based) inside the contig
     * \param ePos      end position (0-based) inside the contig
     * \param rev       \c true if the read is reverse complemented
     */
    Read(const int32_t ctg, const int32_t sPos, const int32_t ePos, const bool rev = false);

    //! Returns contig's identifier.
    /*!
     * \return contig's identifier.
     */
    int32_t getContigId() const;

    //! Returns starting position of the read in the contig
    /*!
     * \return starting position (0-based)
     */
    int32_t getStartPos() const;

    //! Returns ending position of the read in the contig
    /*!
     * \return ending position (0-based)
     */
    int32_t getEndPos() const;

    //! Returns the length of the read.
    /*!
     * \return read's length
     */
    int32_t getLength() const;

    //! Returns whether the read is reverse complemented or not.
    /*!
     * \return \c true if the read is reverse complemented, \c false otherwise.
     */
    bool isReverse();

    bool overlaps( Read &read, int minOverlap = 0 ) const;

    //! Loads a set of (mapped) reads from a bam file.
    /*!
     * Reads with multiple alignments or unmapped are discarded.
     *
     * \param bamReader BamReader object.
     * \param readMap_1 map where the first uniquely mapped pairs are loaded (output)
	 * \param readMap_2 map where the second uniquely mapped pairs are loaded (output)
     * \param coverage vector of coverages of the contigs (output)
	 * \param dupReadMap_1 map where the first mutiple mapped pairs are loaded (output)
	 * \param dupReadMap_2 map where the second mutiple mapped pairs are loaded (output)
	 * \param loadDupReads whether duplicate reads maps should be filled
	 *
     */
    static void loadReadsMap(
        MultiBamReader &bamReader,
        sparse_hash_map< std::string, Read > &readMap_1,
        sparse_hash_map< std::string, Read > &readMap_2,
        std::vector< std::vector<uint32_t> > &coverage,
        bool noMultFilter = false
	);
};

#endif	/* READS_HPP */

