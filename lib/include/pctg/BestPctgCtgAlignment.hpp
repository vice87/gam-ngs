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
    std::vector<MyAlignment> _main;		//!< main alignments (one for each block)
    MyAlignment _left;					//!< left tail alignment
	MyAlignment _right;					//!< right tail alignment

    bool _isCtgReverse; 	//!< whether the contig is reverse complemented or not.
    bool _left_rev;			//!< whether left alignment considers contings in opposite order
    bool _right_rev;		//!< whether right alignment considers contings in opposite order

public:
    //! A constructor.
    /*!
     * \param main an alignment.
     * \param isCtgReverse whether the contig is reverse complemented or not.
     */
	BestPctgCtgAlignment();
    BestPctgCtgAlignment(const MyAlignment &main, const bool isCtgReverse);
	BestPctgCtgAlignment(const std::vector<MyAlignment> &main, const bool isCtgReverse);

	//! A constructor.
	/*!
	 * \param a an alignment.
	 * \param isCtgReverse whether the contig is reverse complemented or not.
	 * \param left left tail alignment
	 * \param right right tail alignment
	 * \param left_rev whether the left alignment has been done considering contigs in opposite order
	 * \param right_rev whether the left alignment has been done considering contigs in opposite order
	 */
	BestPctgCtgAlignment(
		const MyAlignment& main,
		const bool isCtgReverse,
		const MyAlignment& left,
		const MyAlignment& right,
		const bool left_rev = false,
		const bool right_rev = false
	);

	BestPctgCtgAlignment(
		const std::vector<MyAlignment> &main,
		const bool isCtgReverse,
		const MyAlignment &left,
		const MyAlignment &right,
		const bool left_rev = false,
		const bool right_rev = false
	);

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
    const std::vector<MyAlignment>& main() const;
	const MyAlignment& left() const;
	const MyAlignment& right() const;

	double main_homology() const;

	bool is_left_rev() const;
	bool is_right_rev() const;

    //! Tells whether the contig is reverse complemented or not.
    bool isCtgReversed() const;

    //! Assign operator of the BestPctgCtgAlignment class.
    const BestPctgCtgAlignment& operator =(const BestPctgCtgAlignment &orig);

	const MyAlignment& operator[]( const size_t &index ) const;
	uint64_t size() const;
};


#endif	/* BESTPCTGALIGNMENT_HPP */

