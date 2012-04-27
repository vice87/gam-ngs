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
    IdType _aId;        //!< assembly identifier
    IdType _ctgId;      //!< contig identifier
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
     * \param aId       assembly identifier
     * \param ctgId     contig identifier
     * \param size      size of the contig
     * \param position  position of the contig inside a paired contig.
     */
    ContigInPctgInfo(const IdType &aId, const IdType &ctgId, const UIntType &size, const UIntType &position);

    //! A constructor.
    /*!
     * Creates a ContigInPctgInfo object, given a contig ID, its size and a position.
     * \param ctgId     contig identifier
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
     * \param aId assembly identifier
     * \param ctgId contig identifier
     * \param bestAlign a BestPctgCtgAlignment object
     *
     * \sa BestPctgCtgAlignment
     */
    ContigInPctgInfo(const IdType &aId, const IdType &ctgId, const BestPctgCtgAlignment &bestAlign);

    //! A constructor.
    /*!
     * Creates a ContigInPctgInfo object, given a contig ID and a BestPctgCtgAlignment object.
     * \param ctgId contig's identifier
     * \param bestAlign a BestPctgCtgAlignment object
     *
     * \sa BestPctgCtgAlignment
     */
    ContigInPctgInfo(const IdType &ctgId, const BestPctgCtgAlignment &bestAlign);

    //! Gets sequence identifier.
    /*!
     * \return pair of assembly and contig IDs.
     */
    std::pair<IdType,IdType> getId() const;

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

