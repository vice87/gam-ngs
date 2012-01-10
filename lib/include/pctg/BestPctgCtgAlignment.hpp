/*! 
 * \file BestPctgCtgAlignment.hpp
 * \brief Definition of BestPctgCtgAlignment class.
 * \details This file contains the definition of the class representing a best
 * alignment of a contig and a paired contig.
 */

#ifndef BESTPCTGALIGNMENT_HPP
#define	BESTPCTGALIGNMENT_HPP

#include "alignment/my_alignment.hpp"

//! Class implementing a best pctg-contig alignment
class BestPctgCtgAlignment
{

private:
    MyAlignment _a;       //!< an alignment
    bool _isCtgReverse; //!< whether the contig is reverse complemented or not.
    
public:
    //! A constructor.
    /*!
     * \param a an alignment.
     * \param isCtgReverse whether the contig is reverse complemented or not.
     */
    BestPctgCtgAlignment(const MyAlignment& a, const bool isCtgReverse);
    
    //! A copy constructor.
    /*!
     * Creates a copy of a given BestPctgCtgAlignment object.
     * \param orig a best pctg-ctg alignment.
     */
    BestPctgCtgAlignment(const BestPctgCtgAlignment &orig);
    
    //! Gets the alignment.
    /*!
     * \return a reference to the Alignment object.
     */
    const MyAlignment& getAlignment() const;
    
    //! Tells whether the contig is reverse complemented or not.
    const bool& isCtgReversed() const;
    
    //! Assign operator of the BestPctgCtgAlignment class.
    const BestPctgCtgAlignment& operator =(const BestPctgCtgAlignment &orig);
};


#endif	/* BESTPCTGALIGNMENT_HPP */

