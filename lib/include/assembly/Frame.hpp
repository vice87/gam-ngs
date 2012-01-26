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
    int32_t   _ctgId;      //!< contig identifier
    char     _strand;     //!< strand of the frame
    int32_t _begin;      //!< position (0-based) where the block begins in contig
    int32_t _end;        //!< position (0-based) where the block ends in contig
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
    Frame( int32_t ctg, char strand, int32_t begin, int32_t end );
    
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
    
    void setReadsLen( IntType readLen );
    
    void increaseReadsLen( IntType readLen );
    
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
    
    IntType getReadsLen() const;
    
    //! Gets the frame's length.
    /*!
     * \return frame's length.
     */
    int32_t getLength() const;
    
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

