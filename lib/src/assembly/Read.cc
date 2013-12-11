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

#include "OrderingFunctions.hpp"

#include "assembly/Read.hpp"

Read::Read():
        _contigId(0), _startPos(0), _endPos(0), _isRev(false)
{ }

Read::Read(const Read &orig):
        _contigId(orig._contigId), _startPos(orig._startPos),
        _endPos(orig._endPos), _isRev(orig._isRev)
{ }

Read::Read(const int32_t ctg, const int32_t sPos, const int32_t ePos, const bool rev):
        _contigId(ctg), _startPos(sPos), _endPos(ePos), _isRev(rev)
{ }

int32_t Read::getContigId() const
{
    return _contigId;
}

int32_t Read::getStartPos() const
{
    return _startPos;
}

int32_t Read::getEndPos() const
{
    return (_endPos - 1);
}

int32_t Read::getLength() const
{
    return _endPos - _startPos;
}

bool Read::isReverse()
{
    return _isRev;
}

bool Read::overlaps(Read& read, int minOverlap) const
{
    if( this->_contigId != read.getContigId() ) return false;

    int32_t startDiff = this->getEndPos() - read.getStartPos() + 1;
    int32_t endDiff = read.getEndPos() - this->getStartPos() + 1;

    return (startDiff >= minOverlap && endDiff >= minOverlap);
}

void Read::loadReadsMap(
		MultiBamReader &bamReader,
		sparse_hash_map< std::string, Read > &readMap_1,
		sparse_hash_map< std::string, Read > &readMap_2,
		std::vector< std::vector<uint32_t> > &coverage )
{
    // initialize coverage vector
    const RefVector& refVect = bamReader.GetReferenceData();
    coverage.resize( refVect.size() );
    for( uint32_t i=0; i < refVect.size(); i++ ) coverage.at(i).resize( refVect.at(i).RefLength, 0 );

    int32_t nh, xt;
    BamAlignment align;

	bamReader.Rewind();

    while( bamReader.GetNextAlignment(align,true) )
    {
        // discard unmapped reads and reads that have a bad quality
        if( !align.IsMapped() || align.Position < 0 || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

        // se la molteplicità non è stata definita, assumo che sia pari ad 1
        if( !align.GetTag(std::string("NH"),nh) ) nh = 1;
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U';
        bool uniqMapRead = (nh == 1 && xt == 'U');

		if( !uniqMapRead ) continue; // read mappata in modo molteplice

		Read curRead( align.RefID, align.Position, align.GetEndPosition(), align.IsReverseStrand() );

		// insert reads one of the reads hash-tables depending on whether it is the first or second pair
		if( !align.IsPaired() || align.IsFirstMate() ) readMap_1[align.Name] = curRead; else readMap_2[align.Name] = curRead;

		// update vector coverage
		uint32_t read_len = align.GetEndPosition() - align.Position;
		for( int i=0; i < read_len; i++ ) coverage.at(align.RefID).at(align.Position+i) += 1;
    }
}

void Read::loadMasterReadsMap( BamReader &bamMaster, BamReader &bamSlave, sparse_hash_map< std::string, Read > &readMap )
{
    int32_t nh, xt;
    BamAlignment align, masterAlign, prevMasterAlign;

    if( !bamMaster.GetNextAlignment(masterAlign) ) return;

    while( bamSlave.GetNextAlignment(align) )
    {
        if( strnum_cmp( align.Name.c_str(), masterAlign.Name.c_str()) < 0 ) // align.Name < masterAlign.Name
        {
            if( align.Name == prevMasterAlign.Name && align.IsPaired() && align.IsFirstMate() == prevMasterAlign.IsFirstMate() )
            {
                if( !align.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
                if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // l'allineamento non è stato fatto con bwa

                if( nh == 1 && xt == 'U' && align.IsMapped() )
                {
                    if( !prevMasterAlign.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
                    if( !prevMasterAlign.GetTag(std::string("XT"),xt) ) xt = 'U'; // l'allineamento non è stato fatto con bwa
                    if( nh == 1 && xt == 'U' && prevMasterAlign.IsMapped() )
                    {
                        Read masterRead(prevMasterAlign.RefID, prevMasterAlign.Position, prevMasterAlign.GetEndPosition(), prevMasterAlign.IsReverseStrand());
                        std::string readName = prevMasterAlign.Name + (prevMasterAlign.IsFirstMate() ? "1" : "2");
                        readMap.insert( std::make_pair(readName,masterRead) );
                    }
                }
            }

            continue;
        }

        if( strnum_cmp( align.Name.c_str(), masterAlign.Name.c_str()) > 0 ) // align.Name > slaveAlign.Name
        {
            do
            {
                bamMaster.GetNextAlignment(masterAlign);
            }
            while( strnum_cmp(align.Name.c_str(), masterAlign.Name.c_str()) > 0 );
        }

        if( strnum_cmp( align.Name.c_str(), masterAlign.Name.c_str()) == 0 && align.IsPaired() && align.IsFirstMate() != masterAlign.IsFirstMate() )
        {
            prevMasterAlign = masterAlign;
            if( !bamMaster.GetNextAlignment(masterAlign) ) break;
        }

        if( strnum_cmp( align.Name.c_str(), masterAlign.Name.c_str()) == 0 )
        {
            if( !align.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
            if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // l'allineamento non è stato fatto con bwa
            if( nh == 1 && xt == 'U' && align.IsMapped() )
            {
                if( !masterAlign.GetTag(std::string("NH"),nh) ) nh = 1; // se molteplicità in un campo non standard assumo che sia pari ad 1
                if( !masterAlign.GetTag(std::string("XT"),xt) ) xt = 'U'; // l'allineamento non è stato fatto con bwa
                if( nh == 1 && xt == 'U' && masterAlign.IsMapped() )
                {
                    Read masterRead(masterAlign.RefID, masterAlign.Position, masterAlign.GetEndPosition(), masterAlign.IsReverseStrand());
                    std::string readName = masterAlign.Name;

                    if( align.IsPaired() ) readName = readName + (masterAlign.IsFirstMate() ? "1" : "2");

                    readMap.insert( std::make_pair(readName,masterRead) );
                }
            }

            if( !bamMaster.GetNextAlignment(masterAlign) ) break;
        }
    }
}

