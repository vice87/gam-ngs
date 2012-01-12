/*! 
 * \file PairedContig.hpp
 * \brief Definition of PairedContig class.
 * \details This file contains the definition of the class representing a paired 
 * contig, i.e. a contig consisting of one or more contigs merged toghether.
 */

#ifndef PAIREDCONTIG_HPP
#define	PAIREDCONTIG_HPP

#define PCTG_DEFAULT_PREFIX_NAME "PairedContig_"

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <algorithm>

#include "api/BamAux.h"

#include "types.hpp"
#include "assembly/contig.hpp"
#include "pctg/ContigInPctgInfo.hpp"

//! Class implementing a paired contig.
class PairedContig : public Contig
{
private:
    typedef std::map< IdType, ContigInPctgInfo > ContigInfoMap;
    
    ContigInfoMap _masterCtgMap;        //!< map of master contigs in this paired contig
    ContigInfoMap _slaveCtgMap;         //!< map of slave contigs in this paired contig
    IdType _pctgId;                     //!< paired contig ID
    
public:
    //! A constructor.
    /*!
     * This method creates an empty paired contig.
     */
    PairedContig();
    
    //! A constructor.
    /*!
     * Creates a paired contig, given an ID number.
     * \param id identifier of the paired contig
     */
    PairedContig(const IdType &id);
    
    //! A copy constructor.
    /*!
     * Creates copy of a given paired contig.
     * \param orig a paired contig
     */
    PairedContig(const PairedContig &orig);
    
    //! A copy constructor.
    /*!
     * Creates a copy of a given paired contig and a contig. The master/slave
     * map is copied from \c orig, the sequence is copied from \c ctg.
     * The name is reinitilized with the identifier of \c orig.
     * \param orig a paired contig.
     * \param ctg a contig.
     */
    PairedContig(const PairedContig &orig, const Contig &ctg);
    
    
    //! Sets the identifier.
    /*!
     * This method sets the current paired contig identifier.
     * \param id an identifier.
     */
    void setId(const IdType &id);
    
    //! Gets the identifier.
    /*!
     * This method returns the current paired contig identifier.
     * \return the paired contig id.
     */
    IdType getId() const;
    
    //! Gets the master contigs map.
    /*!
     * \return a reference of the master contigs map.
     */
    ContigInfoMap& getMasterCtgMap();
    
    //! Gets the slave contigs map.
    /*!
     * \return a reference of the slave contigs map.
     */
    ContigInfoMap& getSlaveCtgMap();
    
    //! Gets the master contig map.
    /*!
     * \return a constant reference of the master contigs map.
     */
    const ContigInfoMap& getMasterCtgMap() const;
    
    //! Gets the master contig map.
    /*!
     * \return a constant reference of the slave contigs map.
     */
    const ContigInfoMap& getSlaveCtgMap() const;
    
    //! Returns a ContigInPctgInfo object of a contig inside the paired contig.
    /*!
     * \param ctgId a contig identifier.
     * \param isMasterCtg specifies whether \c ctgIs is a master contig or not.
     * \return a ContigInPctgInfo object of the contig.
     * 
     * \sa ContigInPctgInfo
     */
    const ContigInPctgInfo& getContigInfo(const IdType &ctgId, bool isMasterCtg);
    
    //! Gets the starting position of a given ContigInPctgInfo object.
    /*!
     * \param ctgInfo a ContigInPctgInfo object.
     * \return the starting position of \c ctgInfo
     */
    UIntType getContigBegin(const ContigInPctgInfo &ctgInfo) const;
    
    //! Gets the ending position of a given ContigInPctgInfo object.
    /*!
     * \param ctgInfo a ContigInPctgInfo object.
     * \return the ending position of \c ctgInfo
     */
    UIntType getContigEnd(const ContigInPctgInfo &ctgInfo) const;
    
    //! Gets the base position of the contig in the paired contig
    /*!
     * \param ctgId contig's identifier.
     * \param pos offset
     * \param isMasterCtg specifies whether \c ctgId is a master contig or not
     */
    UIntType getBasePosition(
        const IdType &ctgId, 
        const UIntType pos,
        const bool isMasterCtg);
    
    //! Returns whether a ctgId is contained into the master contigs map.
    /*!
     * \param ctgId a contig identifier.
     * \return whether \c ctgId is contained into the master contigs map.
     */
    bool containsMasterCtg(const IdType &ctgId) const;
    
    //! Returns whether a ctgId is contained into the slave contigs map.
    /*!
     * \param ctgId a contig identifier.
     * \return whether \c ctgId is contained into the slave contigs map.
     */
    bool containsSlaveCtg(const IdType &ctgId) const;
    
    //! Assign operator of PairedContig class.
    const PairedContig& operator =(const PairedContig& orig);
};


//! Returns whether a contig of an assembly overlaps a contig of the same assembly in the paired contig.
/*!
 * \param pctg a paired contig
 * \param ctg a contig
 * \param pos position in pctg where the overlap is tested.
 * \param isMasterCtg tells whether \c ctg is a master or slave contig.
 */
bool sameAssemblyCtgsOverlapedBy(const PairedContig &pctg, const Contig &ctg, 
        const UIntType &pos, bool isMasterCtg);


//! Shifts all the contigs inside a paired contig.
/*!
 * \param pctg a paired contig.
 * \param shiftSize shift amount.
 * \return a copy of \c pctg where all the contigs inside it are shifted by \c shiftSize.
 */
PairedContig& shiftOf(PairedContig& pctg, const UIntType& shiftSize);

bool orderPctgsByName(const PairedContig &a, const PairedContig &b);

std::ostream& writePctgDescriptors(std::ostream &os, const std::list<PairedContig> &pctgs, BamTools::RefVector &mcRef, BamTools::RefVector &scRef);

std::ostream& writePctgDescriptor(std::ostream &os, const PairedContig &pctg, BamTools::RefVector &mcRef, BamTools::RefVector &scRef);

#endif	/* PAIREDCONTIG_HPP */

