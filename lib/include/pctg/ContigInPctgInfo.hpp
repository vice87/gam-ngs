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

//! Class storing the properties of a contig inside a paired contig.
class ContigInPctgInfo
{
    
private:
    IdType _ctgId;      //!< contig's identifier
    UIntType _size;     //!< size of the contig
    UIntType _position; //!< position inside the paired contig
    UIntType _gaps;     //!< number of gaps
    bool _reversed;     //!< whether the contig is reverse complemented or not.
    
public:
    //! A constructor.
    /*!
     * Creates an empty ContigInPctgInfo object.
     */
    ContigInPctgInfo();
    
    //! A constructor.
    /*!
     * Creates a ContigInPctgInfo object, given a contig ID, its size and a position.
     * \param ctgId     contig's identifier
     * \param size      size of the contig
     * \param position  position of the contig inside a paired contig.
     */
    ContigInPctgInfo(const IdType &ctgId, const UIntType &size, const UIntType &position);
    
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
    ContigInPctgInfo(const IdType &ctgId, const BestPctgCtgAlignment &bestAlign);
    
    //! Gets contig's identifier.
    /*!
     * \return the contig's ID.
     */
    IdType getId() const;
    
    //! Gets the position of the first nucleotide of the contig inside the paired contig.
    /*!
     * \return position of the first nucleotide.
     */
    UIntType getFirstNucleotidePos() const;
    
    //! Gets the position of the last nucleotide of the contig inside the paired contig.
    /*!
     * \return position of the last nucleotide.
     */
    UIntType getLastNucleotidePos() const;
    
    //! Tells whether the contig is reverse complemented or not.
    /*!
     * \return whether the contig is reverse complemented or not.
     */
    const bool& isReversed() const;
    
    //! Sets the amount of gaps.
    /*!
     * \param gaps number of gaps.
     */
    void setGaps(const UIntType& gaps); 
    
    //! Sets the position.
    /*!
     * \param pos position of the contig inside a paired contig.
     */
    void setPosition(const UIntType& pos);
};


#endif	/* CONTIGINPCTGINFO_HPP */

