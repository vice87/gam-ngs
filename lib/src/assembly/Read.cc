/* 
 * File:   Read.cc
 * Author: Riccardo Vicedomini
 * 
 * Created on 22 maggio 2011, 14.33
 */

#include "OrderingFunctions.hpp"

#include "assembly/Read.hpp"

Read::Read():
        _contigId(0), _startPos(0), _endPos(0), _isRev(false)
{ }

Read::Read(const Read &orig):
        _contigId(orig._contigId), _startPos(orig._startPos),
        _endPos(orig._endPos), _isRev(orig._isRev)
{ }

Read::Read(const IdType ctg, const IdType sPos, const IdType ePos, const bool rev):
        _contigId(ctg), _startPos(sPos), _endPos(ePos), _isRev(rev)
{ }

IdType Read::getContigId() const
{
    return _contigId;
}

IdType Read::getStartPos() const
{
    return _startPos;
}

IdType Read::getEndPos() const
{
    return (_endPos - 1);
}

IntType Read::getLength() const
{
    return _endPos - _startPos;
}

bool Read::isReverse()
{
    return _isRev;
}

void Read::getReadMap( BamReader &bamMaster, BamReader &bamSlave, sparse_hash_map< std::string, Read > &readMap )
{    
    uint32_t nh;
    BamAlignment align, slaveAlign, prevSlaveAlign;
        
    readMap.clear();
    prevSlaveAlign.Name = "";
    
    if( !bamSlave.GetNextAlignment(slaveAlign) ) return;
    
    while( bamMaster.GetNextAlignment(align) )
    {
        if( strnum_cmp( align.Name.c_str(), slaveAlign.Name.c_str()) < 0 ) // align.Name < slaveAlign.Name
        {            
            if( align.Name == prevSlaveAlign.Name && align.IsPaired() && align.IsFirstMate() == prevSlaveAlign.IsFirstMate() )
            {
                if( !align.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
                if( nh == 1 && align.IsMapped() )
                {
                    if( !prevSlaveAlign.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
                    if( nh == 1 && prevSlaveAlign.IsMapped() )
                    {
                        Read slaveRead(IdType(prevSlaveAlign.RefID), IdType(prevSlaveAlign.Position), IdType(prevSlaveAlign.GetEndPosition()), prevSlaveAlign.IsReverseStrand());
                        std::string readName = prevSlaveAlign.Name + (prevSlaveAlign.IsFirstMate() ? "1" : "2");
                        readMap.insert( std::make_pair(readName,slaveRead) );
                    }
                }
            }
            
            continue;
        }
        
        if( strnum_cmp( align.Name.c_str(), slaveAlign.Name.c_str()) > 0 ) // align.Name > slaveAlign.Name
        {
            do
            {
                bamSlave.GetNextAlignment(slaveAlign);
            }
            while( strnum_cmp(align.Name.c_str(), slaveAlign.Name.c_str()) > 0 );
        }
        
        if( strnum_cmp( align.Name.c_str(), slaveAlign.Name.c_str()) == 0 && align.IsPaired() && align.IsFirstMate() != slaveAlign.IsFirstMate() )
        {            
            prevSlaveAlign = slaveAlign;
            if( !bamSlave.GetNextAlignment(slaveAlign) ) break;
        }
        
        if( strnum_cmp( align.Name.c_str(), slaveAlign.Name.c_str()) == 0 )
        {
            if( !align.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
            if( nh == 1 && align.IsMapped() )
            {
                if( !slaveAlign.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
                if( nh == 1 && slaveAlign.IsMapped() )
                {
                    Read slaveRead(IdType(slaveAlign.RefID), IdType(slaveAlign.Position), IdType(slaveAlign.GetEndPosition()), slaveAlign.IsReverseStrand());
                    std::string readName = slaveAlign.Name;
                    
                    if( align.IsPaired() )
                    {
                        readName = readName + (slaveAlign.IsFirstMate() ? "1" : "2");
                    }

                    readMap.insert( std::make_pair(readName,slaveRead) );
                }
            }
            
            if( !bamSlave.GetNextAlignment(slaveAlign) ) break;
        }
    }
}

