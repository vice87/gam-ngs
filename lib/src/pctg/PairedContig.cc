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
        _masterCtgMap(orig._masterCtgMap), _slaveCtgMap(orig._slaveCtgMap)
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


bool PairedContig::containsMasterCtg(const IdType& ctgId) const
{
    return (this->_masterCtgMap).find( ctgId ) != (this->_masterCtgMap).end();
}


bool PairedContig::containsSlaveCtg(const IdType& ctgId) const
{
    return (this->_slaveCtgMap).find( ctgId ) != (this->_slaveCtgMap).end();
}

const ContigInPctgInfo& PairedContig::getContigInfo(const IdType &ctgId, bool isMasterCtg)
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


UIntType PairedContig::getBasePosition(
        const IdType &ctgId, 
        const UIntType pos,
        const bool isMasterCtg)
{
    ContigInPctgInfo ctgInfo( this->getContigInfo(ctgId,isMasterCtg) );
    
    if( !ctgInfo.isReversed() ) return ctgInfo.getFirstNucleotidePos() + pos;
    
    return ctgInfo.getLastNucleotidePos() - pos;
}


const PairedContig& PairedContig::operator =(const PairedContig& orig)
{
    ((Contig *)this)->operator =((Contig)orig);
    this->_pctgId = orig._pctgId;
    this->_masterCtgMap = orig._masterCtgMap;
    this->_slaveCtgMap = orig._slaveCtgMap;
    
    return *this;
}


bool sameAssemblyCtgsOverlapedBy(const PairedContig &pctg, const Contig &ctg, 
        const UIntType &pos, bool isMasterCtg)
{
	return false;
	if( !isMasterCtg ) return false;

    typedef std::map< IdType, ContigInPctgInfo > ContigInfoMap;
    
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
    
    return false;
}


PairedContig& shiftOf(PairedContig& pctg, const UIntType& shiftSize)
{
    typedef std::map< IdType, ContigInPctgInfo > ContigInfoMap;
    
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


std::ostream& writePctgDescriptors( std::ostream &os, const std::list<PairedContig> &pctgs, BamTools::RefVector &mcRef, BamTools::RefVector &scRef )
{
    os << "#Name\tSize\tAssembly\tContigID\tBeginInPctg\tEndInPctg" << std::endl;
    
    std::list< PairedContig >::const_iterator i;
    for( i = pctgs.begin(); i != pctgs.end(); i++ ) writePctgDescriptor(os,*i, mcRef, scRef);
    
    return os;
}


std::ostream& writePctgDescriptor( std::ostream &os, const PairedContig &pctg, BamTools::RefVector &mcRef, BamTools::RefVector &scRef )
{
    typedef std::map< IdType, ContigInPctgInfo > ContigInfoMap;
    
    ContigInfoMap masterCtgs = pctg.getMasterCtgMap();
    ContigInfoMap slaveCtgs = pctg.getSlaveCtgMap();
    
    ContigInfoMap::const_iterator ctg;
    for( ctg = masterCtgs.begin(); ctg != masterCtgs.end(); ctg++ )
    {
        os << pctg.name() << "\t"
                << pctg.size() << "\t"
                << "Master" << "\t"
                << mcRef.at( ctg->first ).RefName << "\t"
                << pctg.getContigBegin( ctg->second ) << "\t"
                << pctg.getContigEnd( ctg->second ) << "\t"
                << std::endl;
    }
    
    for( ctg = slaveCtgs.begin(); ctg != slaveCtgs.end(); ctg++ )
    {
        os << pctg.name() << "\t"
                << pctg.size() << "\t"
                << "Slave" << "\t"
                << scRef.at( ctg->first ).RefName << "\t"
                << pctg.getContigBegin( ctg->second ) << "\t"
                << pctg.getContigEnd( ctg->second ) << "\t"
                << std::endl;
    }
    
    return os;
}
