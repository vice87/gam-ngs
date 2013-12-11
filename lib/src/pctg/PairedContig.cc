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

#include "pctg/PairedContig.hpp"

PairedContig::PairedContig() : _pctgId(0)
{}


PairedContig::PairedContig(const IdType& id) : _pctgId(id)
{
    std::ostringstream ss;
    ss << id;
    this->set_name( std::string(PCTG_DEFAULT_PREFIX_NAME) + ss.str() );
}


PairedContig::PairedContig(const PairedContig& orig):
        Contig((Contig)orig), _pctgId(orig._pctgId),
        _masterCtgMap(orig._masterCtgMap), _slaveCtgMap(orig._slaveCtgMap),
        _masterCtgs(orig._masterCtgs), _slaveCtgs(orig._slaveCtgs),
        _mergeList(orig._mergeList), _dupRegionsEst(orig._dupRegionsEst)
{
    std::ostringstream ss;
    ss << orig._pctgId;
    this->set_name( std::string(PCTG_DEFAULT_PREFIX_NAME) + ss.str() );
}


PairedContig::PairedContig(const PairedContig& orig, const Contig& ctg):
        Contig(ctg)
{
    this->_pctgId = orig._pctgId;
    this->_masterCtgMap = orig._masterCtgMap;
    this->_slaveCtgMap = orig._slaveCtgMap;
	this->_mergeList = orig._mergeList;

	this->_masterCtgs = orig._masterCtgs;
	this->_slaveCtgs = orig._slaveCtgs;

	this->_dupRegionsEst = orig._dupRegionsEst;
    this->set_name(orig.name());
}


void PairedContig::setId(const IdType &id)
{
    this->_pctgId = id;

    std::ostringstream ss;
    ss << id;
    this->set_name( std::string(PCTG_DEFAULT_PREFIX_NAME) + ss.str() );
}

IdType PairedContig::getId() const
{
    return this->_pctgId;
}

PairedContig::ContigInfoMap& PairedContig::getMasterCtgMap()
{
    return this->_masterCtgMap;
}

PairedContig::ContigInfoMap& PairedContig::getSlaveCtgMap()
{
    return this->_slaveCtgMap;
}

const PairedContig::ContigInfoMap& PairedContig::getMasterCtgMap() const
{
    return this->_masterCtgMap;
}

const PairedContig::ContigInfoMap& PairedContig::getSlaveCtgMap() const
{
    return this->_slaveCtgMap;
}


std::list< CtgInPctgInfo >& PairedContig::getMergeList()
{
	return _mergeList;
}

const std::list< CtgInPctgInfo >& PairedContig::getMergeList() const
{
	return _mergeList;
}


bool PairedContig::containsMasterCtg(const int32_t ctgId) const
{
    return (this->_masterCtgMap).find(ctgId) != (this->_masterCtgMap).end();
}


bool PairedContig::containsSlaveCtg(const int32_t ctgId) const
{
    return (this->_slaveCtgMap).find(ctgId) != (this->_slaveCtgMap).end();
}

const ContigInPctgInfo& PairedContig::getContigInfo(const int32_t ctgId, bool isMasterCtg) const
{
    if(isMasterCtg) return this->_masterCtgMap.find(ctgId)->second;
    return this->_slaveCtgMap.find(ctgId)->second;
}

ContigInPctgInfo& PairedContig::getContigInfo(const int32_t ctgId, bool isMasterCtg)
{
	if(isMasterCtg) return this->_masterCtgMap.find(ctgId)->second;
	return this->_slaveCtgMap.find(ctgId)->second;
}


UIntType PairedContig::getContigBegin(const ContigInPctgInfo& ctgInfo) const
{
    if( ctgInfo.isReversed() ) return ctgInfo.getLastNucleotidePos();

    return ctgInfo.getFirstNucleotidePos();
}


UIntType PairedContig::getContigEnd(const ContigInPctgInfo& ctgInfo) const
{
    if( ctgInfo.isReversed() ) return ctgInfo.getFirstNucleotidePos();

    return ctgInfo.getLastNucleotidePos();
}


uint64_t PairedContig::getBasePosition(
        const int32_t ctgId,
        const UIntType pos,
        const bool isMasterCtg)
{
    const ContigInPctgInfo &ctgInfo = this->getContigInfo(ctgId,isMasterCtg);

	int64_t bp = ctgInfo.getFirstNucleotidePos() + int64_t(pos);
	if( !ctgInfo.isReversed() ) return (bp < 0 ? 0 : uint64_t(bp));

    //return ctgInfo.getLastNucleotidePos() - pos;
    int64_t lnp = ctgInfo.getLastNucleotidePos();
    return ( lnp >= pos ? lnp-pos : 0 );
}


const std::set<int32_t>& PairedContig::getMasterCtgIdSet()
{
	return (this->_masterCtgs);
}


void PairedContig::addMasterCtgId( int32_t id )
{
	(this->_masterCtgs).insert(id);
}

void PairedContig::addSlaveCtgId( int32_t id )
{
	(this->_slaveCtgs).insert(id);
}


void PairedContig::addDupRegion( uint64_t length, uint64_t ctgId )
{
	(this->_dupRegionsEst).push_back( std::make_pair(length,ctgId) );
}


const std::list< std::pair<uint64_t,uint64_t> >& PairedContig::getDupRegions()
{
	return this->_dupRegionsEst;
}


const PairedContig& PairedContig::operator =(const PairedContig& orig)
{
    ((Contig *)this)->operator =((Contig)orig);
    this->_pctgId = orig._pctgId;
    this->_masterCtgMap = orig._masterCtgMap;
    this->_slaveCtgMap = orig._slaveCtgMap;

	this->_masterCtgs = orig._masterCtgs;
	this->_slaveCtgs = orig._slaveCtgs;

	this->_mergeList = orig._mergeList;
	this->_dupRegionsEst = orig._dupRegionsEst;

    return *this;
}


bool sameAssemblyCtgsOverlapedBy(const PairedContig &pctg, const Contig &ctg,
        const UIntType &pos, bool isMasterCtg)
{
	return false;
	if( !isMasterCtg ) return false;

   /* typedef std::map< IdType, ContigInPctgInfo > ContigInfoMap;

    UIntType cBegin = pos, cEnd = pos + ctg.size();

    ContigInfoMap ctgMap = (isMasterCtg) ? pctg.getMasterCtgMap() : pctg.getSlaveCtgMap();

    ContigInfoMap::const_iterator i;
    for(i = ctgMap.begin(); i != ctgMap.end(); i++ )
    {
        UIntType iBegin = pctg.getContigBegin(i->second);
        UIntType iEnd = pctg.getContigEnd(i->second);

        if( iBegin > iEnd ) std::swap(iBegin,iEnd);

        if( iEnd >= cBegin && cEnd >= iBegin)
	{
		if( iBegin >= cBegin )
		{
			if( iEnd >= cEnd )
			{
				if( 100*(cEnd-iBegin+1) > 25*(iEnd-iBegin+1) || 100*(cEnd-iBegin+1) > 25*(cEnd-cBegin+1))
				return true;
			}
			else
			{
				if( 100*(iEnd-iBegin+1) > 25*(iEnd-iBegin+1) || 100*(iEnd-iBegin+1) > 25*(cEnd-cBegin+1))
				return true;
			}
		} else {
			if( iEnd >= cEnd )
			{
				if( 100*(cEnd-cBegin+1) > 25*(iEnd-iBegin+1) || 100*(cEnd-cBegin+1) > 25*(cEnd-cBegin+1))
                                return true;
			}
			else
			{
				if( 100*(iEnd-cBegin+1) > 25*(iEnd-iBegin+1) || 100*(iEnd-cBegin+1) > 25*(cEnd-cBegin+1))
                                return true;
			}

		}
	}
    }

    return false; **/
}


PairedContig& shiftOf(PairedContig& pctg, const UIntType& shiftSize)
{
    typedef std::map< int32_t, ContigInPctgInfo > ContigInfoMap;

    //PairedContig out(pctg);
    ContigInfoMap::iterator j;

    for(j = pctg.getMasterCtgMap().begin(); j != pctg.getMasterCtgMap().end(); j++)
    {
        (j->second).setPosition((j->second).getFirstNucleotidePos() + shiftSize);
    }

    for(j = pctg.getSlaveCtgMap().begin(); j != pctg.getSlaveCtgMap().end(); j++)
    {
        (j->second).setPosition((j->second).getFirstNucleotidePos() + shiftSize);
    }

    return pctg;
}


bool orderPctgsByName(const PairedContig &a, const PairedContig &b)
{
    int num_a = atoi( a.name().substr( std::string(PCTG_DEFAULT_PREFIX_NAME).size() ).c_str() );
    int num_b = atoi( b.name().substr( std::string(PCTG_DEFAULT_PREFIX_NAME).size() ).c_str() );

    return num_a < num_b;
}


std::ostream& writePctgDescriptors(
	std::ostream &os,
	const std::list<PairedContig> &pctgs,
	RefSequence &masterRef,
	RefSequence &slaveRef,
	uint64_t pctg_id )
{
    os << "#Name\tSize\tAssembly\tContigID\tBegin\tEnd\tReversed" << std::endl;

	uint64_t j=0;

    std::list< PairedContig >::const_iterator i;
    for( i = pctgs.begin(); i != pctgs.end(); i++ )
	{
		if(j==pctg_id) os << "# ----------------------------------------------------" << std::endl;

		writePctgDescriptor(os,*i, masterRef, slaveRef);
		j++;
	}

    return os;
}


std::ostream& writePctgDescriptor(
	std::ostream &os,
	const PairedContig &pctg,
	RefSequence &masterRef,
	RefSequence &slaveRef )
{
	const std::list< CtgInPctgInfo >& mergeList = pctg.getMergeList();

	for( std::list<CtgInPctgInfo>::const_iterator it = mergeList.begin(); it != mergeList.end(); it++ )
	{
		os << pctg.name() << "\t"
		<< pctg.size() << "\t"
		<< (it->isMaster() ? "Master" : "Slave") << "\t"
		<< (it->isMaster() ? masterRef[it->getId()].RefName : slaveRef[it->getId()].RefName)  << "\t"
		<< it->getStart() << "\t"
		<< it->getEnd() << "\t"
		<< (it->isReversed() ? "R" : "F")
		<< std::endl;
	}

	return os;


    typedef std::map< int32_t, ContigInPctgInfo > ContigInfoMap;

    ContigInfoMap masterCtgs = pctg.getMasterCtgMap();
    ContigInfoMap slaveCtgs = pctg.getSlaveCtgMap();

    ContigInfoMap::const_iterator ctg;
    for( ctg = masterCtgs.begin(); ctg != masterCtgs.end(); ctg++ )
    {
        os << pctg.name() << "\t"
                << pctg.size() << "\t"
                << "Master" << "\t"
                << masterRef[ctg->first].RefName << "\t"
                << pctg.getContigBegin( ctg->second ) << "\t"
                << pctg.getContigEnd( ctg->second ) << "\t"
                << std::endl;
    }

    for( ctg = slaveCtgs.begin(); ctg != slaveCtgs.end(); ctg++ )
    {
        os << pctg.name() << "\t"
                << pctg.size() << "\t"
                << "Slave" << "\t"
                << slaveRef[ctg->first].RefName << "\t"
                << pctg.getContigBegin( ctg->second ) << "\t"
                << pctg.getContigEnd( ctg->second ) << "\t"
                << std::endl;
    }

    return os;
}
