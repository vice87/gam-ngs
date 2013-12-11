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
 * \file Frame.hpp
 * \brief Definition of Frame class.
 * \details This file contains the definition of the class representing a frame
 *          in an assembly.
 */

#ifndef FRAME_HPP
#define	FRAME_HPP

#include "types.hpp"
#include "assembly/Read.hpp"

//! Class implementing a frame.
class Frame
{
    //! Inserts a frame into an output stream.
    friend std::ostream &operator<<( std::ostream &output, const Frame &frame );

    //! Extracts a frame from an input stream.
    friend std::istream &operator>>( std::istream &input, Frame &frame );

private:
    int32_t _assemblyId;        //!< assembly identifier
    int32_t _ctgId;             //!< contig identifier
    char    _strand;            //!< strand of the frame
    int32_t _begin;             //!< position (0-based) where the block begins in contig
    int32_t _end;               //!< position (0-based) where the block ends in contig
    UIntType _readsLen;         //!< sum of the lengths of all reads overlapping the frame.
    UIntType _blockReadsLen;    //!< sum of the lengths of the reads belonging to the block in the frame.

public:
    //! A constructor with no arguments.
    Frame();

    //! A copy constructor.
    /*!
     * Creates a copy of a frame
     * \param orig a Frame object
     */
    Frame( const Frame& orig );

    //! A constructor which initializes frame's attributes.
    /*!
     * \param ctg       contig identifier
     * \param strand    strand of the frame (\c '+' or \c '-')
     * \param begin     beginning of the frame in contig (0-based)
     * \param length    length of the frame
     */
    Frame( int32_t ctg, char strand, int32_t begin, int32_t end );

    //! Sets frame's assembly ID
    /*!
     * \param id    assembly identifier
     */
    void setAssemblyId( int32_t aId );

    //! Sets frame's strand
    /*!
     * \param strand    strand of the frame.
     */
    void setStrand( char strand );

    //! Sets the beginning of the frame.
    /*!
     * \param begin position where the frame starts.
     */
    void setBegin( int32_t begin );

    //! Sets the end of the frame.
    /*!
     * \param end position where the frame ends.
     */
    void setEnd( int32_t end );

    void setReadsLen( UIntType readLen );

    void increaseReadsLen( UIntType readLen );

    void setBlockReadsLen( UIntType len );

    void increaseBlockReadsLen( UIntType len );

    //! Gets assembly identifier.
    /*!
     * \return the identifier of the assembly.
     */
    int32_t getAssemblyId() const;

    //! Gets contig identifier.
    /*!
     * \return the identifier of the contig containing the frame.
     */
    int32_t getContigId() const;

    //! Gets frame's strand.
    /*!
     * \return strand of the frame.
     */
    char getStrand() const;

    //! Gets the position (0-based) where the frame begins in contig
    /*!
     * \return the position where the frame starts.
     */
    int32_t getBegin() const;

    //! Gets the position (0-based) where the frame ends in contig
    /*!
     * \return the position where the frame ends.
     */
    int32_t getEnd() const;

    UIntType getReadsLen() const;

    UIntType getBlockReadsLen() const;

    //! Gets the frame's length.
    /*!
     * \return frame's length.
     */
    int32_t getLength() const;

    //! Returns whether a read overlaps the frame
    /*!
     * \param read a read
     * \param minOverlap number of bases required for a read to overlap the frame
     *
     * \return true if the read overlap the frame, false otherwise
     */
    bool overlaps( Read &read, int minOverlap = 0 ) const;

    //! Overloaded assig operator for the Frame class.
    const Frame& operator=(const Frame& frame);

    //! Less-than operator for the Frame class.
    bool operator <(const Frame &frame) const;

    //! Equality operator for the Frame class.
    bool operator ==( const Frame &frame ) const;

    //! Returns whether two frames overlap.
    static bool frameOverlap( const Frame &a, const Frame &b, double minOverlap = 85.0 );
};

#endif	/* FRAME_HPP */

