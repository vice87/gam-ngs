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
#include <set>
#include <algorithm>

#include "api/BamAux.h"

#include "types.hpp"
#include "assembly/contig.hpp"
#include "assembly/RefSequence.hpp"
#include "pctg/ContigInPctgInfo.hpp"
#include "pctg/CtgInPctgInfo.hpp"

//! Class implementing a paired contig.
class PairedContig : public Contig
{
private:
    typedef std::map< int32_t, ContigInPctgInfo > ContigInfoMap;

    ContigInfoMap _masterCtgMap;        //!< map of master contigs in this paired contig
    ContigInfoMap _slaveCtgMap;         //!< map of slave contigs in this paired contig
    IdType _pctgId;                     //!< paired contig ID

    std::list< CtgInPctgInfo > _mergeList;
	std::set< int32_t > _masterCtgs;
	std::set< int32_t > _slaveCtgs;

    std::list< std::pair<uint64_t,uint64_t> > _dupRegionsEst; //!< list of possible duplicate regions lengths

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
    std::set< int32_t >& getMasterIds() { return _masterCtgs; }

    //! Gets the slave contigs map.
    /*!
     * \return a reference of the slave contigs map.
     */
    ContigInfoMap& getSlaveCtgMap();
    std::set< int32_t >& getSlaveIds() { return _slaveCtgs; }

    //! Gets the master contig map.
    /*!
     * \return a constant reference of the master contigs map.
     */
    const ContigInfoMap& getMasterCtgMap() const;
    const std::set< int32_t >& getMasterIds() const { return _masterCtgs; }

    //! Gets the master contig map.
    /*!
     * \return a constant reference of the slave contigs map.
     */
    const ContigInfoMap& getSlaveCtgMap() const;
    const std::set< int32_t >& getSlaveIds() const { return _slaveCtgs; }

	std::list< CtgInPctgInfo >& getMergeList();
	const std::list< CtgInPctgInfo >& getMergeList() const;

    //! Returns a ContigInPctgInfo object of a contig inside the paired contig.
    /*!
     * \param ctgId a contig identifier.
     * \param isMasterCtg specifies whether \c ctgIs is a master contig or not.
     * \return a ContigInPctgInfo object of the contig.
     *
     * \sa ContigInPctgInfo
     */
    const ContigInPctgInfo& getContigInfo(const int32_t ctgId, bool isMasterCtg) const;
	ContigInPctgInfo& getContigInfo(const int32_t ctgId, bool isMasterCtg);

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
    uint64_t getBasePosition(
        const int32_t ctgId,
        const UIntType pos,
        const bool isMasterCtg);

	void addDupRegion( uint64_t length, uint64_t ctgId );

	void addMasterCtgId( int32_t id );
	void addSlaveCtgId( int32_t id );

	const std::set<int32_t>& getMasterCtgIdSet();

	const std::list< std::pair<uint64_t,uint64_t> >& getDupRegions();

    //! Returns whether a ctgId is contained into the master contigs map.
    /*!
     * \param ctgId a contig identifier.
     * \return whether \c ctgId is contained into the master contigs map.
     */
    bool containsMasterCtg(const int32_t ctgId) const;

    //! Returns whether a ctgId is contained into the slave contigs map.
    /*!
     * \param ctgId a contig identifier.
     * \return whether \c ctgId is contained into the slave contigs map.
     */
    bool containsSlaveCtg(const int32_t ctgId) const;

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

std::ostream& writePctgDescriptors
(
	std::ostream &os,
	const std::list<PairedContig> &pctgs,
	RefSequence &masterRef,
	RefSequence &slaveRef,
	uint64_t pctg_id
);

std::ostream& writePctgDescriptor
(
	std::ostream &os,
	const PairedContig &pctg,
	RefSequence &masterRef,
	RefSequence &slaveRef
);

#endif	/* PAIREDCONTIG_HPP */

