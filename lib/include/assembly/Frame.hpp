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
    IdType   _ctgId;      //!< contig identifier
    char     _strand;     //!< strand of the frame
    UIntType _begin;      //!< position (0-based) where the block begins in contig
    UIntType _end;        //!< position (0-based) where the block ends in contig
    IntType _readsLen;   //!< sum of the lengths of block's reads

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
    Frame( IdType ctg, char strand, UIntType begin, UIntType end );
    
    //! Sets frame's strand
    /*!
     * \param strand    strand of the frame.
     */    
    void setStrand( char strand );
    
    //! Sets the beginning of the frame.
    /*!
     * \param begin position where the frame starts.
     */
    void setBegin( UIntType begin );
    
    //! Sets the end of the frame.
    /*!
     * \param end position where the frame ends.
     */
    void setEnd( UIntType end );
    
    void setReadsLen( IntType readLen );
    
    void increaseReadsLen( IntType readLen );
    
    //! Gets contig identifier.
    /*!
     * \return the identifier of the contig containing the frame.
     */
    IdType getContigId() const;
    
    //! Gets frame's strand.
    /*!
     * \return strand of the frame.
     */
    char getStrand() const;
    
    //! Gets the position (0-based) where the frame begins in contig
    /*!
     * \return the position where the frame starts.
     */
    UIntType getBegin() const;
    
    //! Gets the position (0-based) where the frame ends in contig
    /*!
     * \return the position where the frame ends.
     */
    UIntType getEnd() const;
    
    IntType getReadsLen() const;
    
    //! Gets the frame's length.
    /*!
     * \return frame's length.
     */
    UIntType getLength() const;
    
    //! Overloaded assig operator for the Frame class.
    const Frame& operator=(const Frame& frame); 
    
    //! Less-than operator for the Frame class.
    bool operator <(const Frame &frame) const;
    
    //! Equality operator for the Frame class.
    bool operator ==( const Frame &frame ) const;
    
    //! Returns whether two frames overlap.
    static bool frameOverlap( const Frame &a, const Frame &b, double minOverlap = 50.0 );
};

#endif	/* FRAME_HPP */

