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

#include <cstdlib>
#include <iostream>
#include "assembly/Frame.hpp"

Frame::Frame()
        : _assemblyId(0), _ctgId(0), _strand('?'), _begin(0), _end(0), _readsLen(0), _blockReadsLen(0)
{}

Frame::Frame( int32_t ctg, char strand, int32_t begin, int32_t end )
        : _assemblyId(0), _ctgId(ctg), _strand(strand), _begin(begin), _end(end), _readsLen(0), _blockReadsLen(0)
{}

Frame::Frame(const Frame& orig)
        : _assemblyId(orig._assemblyId), _ctgId(orig._ctgId), _strand(orig._strand),
          _begin(orig._begin),_end(orig._end), _readsLen(orig._readsLen), _blockReadsLen(orig._blockReadsLen)
{}


void Frame::setAssemblyId( int32_t aId )
{
    _assemblyId = aId;
}

void Frame::setStrand( char strand )
{
    _strand = strand;
}

void Frame::setBegin( int32_t begin )
{
    _begin = begin;
}

void Frame::setEnd( int32_t end )
{
    _end = end;
}

void Frame::setReadsLen( UIntType readLen )
{
    _readsLen = readLen;
}

void Frame::increaseReadsLen( UIntType readLen )
{
    _readsLen += readLen;
}

void Frame::setBlockReadsLen(UIntType len)
{
    _blockReadsLen = len;
}

void Frame::increaseBlockReadsLen(UIntType len)
{
    _blockReadsLen += len;
}

int32_t Frame::getAssemblyId() const
{
    return _assemblyId;
}

int32_t Frame::getContigId() const
{
    return _ctgId;
}

char Frame::getStrand() const
{
    return _strand;
}

int32_t Frame::getBegin() const
{
    return _begin;
}

int32_t Frame::getEnd() const
{
    return _end;
}

UIntType Frame::getReadsLen() const
{
    return _readsLen;
}

UIntType Frame::getBlockReadsLen() const
{
    return _blockReadsLen;
}

int32_t Frame::getLength() const
{
    return (_end < _begin) ? 0 : _end - _begin + 1;
}


bool Frame::overlaps(Read& read, int minOverlap) const
{
    if( this->getContigId() != read.getContigId() ) return false;

	IntType startDiff = IntType(this->getEnd()) - IntType(read.getStartPos()) + 1;
    IntType endDiff = IntType(read.getEndPos()) - IntType(this->getBegin()) + 1;

    return (startDiff >= minOverlap && endDiff >= minOverlap);
}


bool Frame::frameOverlap( const Frame& a, const Frame& b, double minOverlap )
{
    if( a.getContigId() != b.getContigId() ) return false;

	int32_t end = std::min( a.getEnd(), b.getEnd() );

    if( a.getBegin() <= b.getBegin() )
    {
        double overlap = (a.getEnd()>=b.getBegin()) ? (100.0 * (double)(end - b.getBegin() + 1) / b.getLength()) : 0.0;
        return ( overlap >= minOverlap );
    }
    else
    {
        double overlap = (b.getEnd()>=a.getBegin()) ? (100.0 * (double)(end - a.getBegin() + 1) / b.getLength()) : 0.0;
        return ( overlap >= minOverlap );
    }
}


const Frame&
Frame::operator=(const Frame& frame)
{
    this->_assemblyId = frame._assemblyId;
    this->_ctgId = frame._ctgId;
    this->_begin = frame._begin;
    this->_end = frame._end;
    this->_strand = frame._strand;
    this->_readsLen = frame._readsLen;
    this->_blockReadsLen = frame._blockReadsLen;

    return *this;
}


bool
Frame::operator <(const Frame& frame) const
{
    return( this->_assemblyId < frame.getAssemblyId() ||
            (this->_assemblyId == frame.getAssemblyId() && this->_ctgId < frame.getContigId()) ||
            (this->_assemblyId == frame.getAssemblyId() && this->_ctgId == frame.getContigId() && this->_begin < frame.getBegin()) ||
            (this->_assemblyId == frame.getAssemblyId() && this->_ctgId == frame.getContigId() && this->_begin == frame.getBegin() && this->_end < frame.getEnd())
          );
}


bool
Frame::operator ==(const Frame& frame) const
{
    return( this->_assemblyId == frame.getAssemblyId() &&
            this->_ctgId == frame.getContigId() &&
            this->_strand == frame.getStrand() &&
            this->_begin == frame.getBegin() &&
            this->_end == frame.getEnd() );
}


std::ostream &operator<<( std::ostream &output, const Frame &frame )
{
    output << 0 << "\t"
           << frame.getContigId() << "\t" << frame.getStrand() << "\t"
           << frame.getBegin() << "\t" << frame.getEnd() << "\t"
           << frame.getBlockReadsLen() << "\t" << frame.getReadsLen();

    return output;
}

std::istream &operator>>( std::istream &input, Frame &frame )
{
    int32_t assemblyId; // unused
    int32_t ctgId;
    char strand;
    int32_t begin, end;
    UIntType readsLen, blockReadsLen;

    input >> assemblyId >> ctgId >> strand >> begin >> end >> blockReadsLen >> readsLen;

    frame = Frame( ctgId, strand, begin, end );
    frame.setReadsLen(readsLen);
    frame.setBlockReadsLen(blockReadsLen);

    return input;
}