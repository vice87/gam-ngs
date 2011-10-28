/* 
 * File:   Frame.cc
 * Author: Riccardo Vicedomini
 * 
 * Created on 22 maggio 2011, 17.40
 */

#include <cstdlib>
#include <iostream>
#include "assembly/Frame.hpp"

Frame::Frame() {}

Frame::Frame( IdType ctg, char strand, UIntType begin, UIntType end )
        : _ctgId(ctg), _strand(strand), _begin(begin), _end(end) {}

Frame::Frame(const Frame& orig)
        : _ctgId(orig._ctgId), _strand(orig._strand), _begin(orig._begin),
          _end(orig._end) {}


void Frame::setStrand( char strand )
{
    _strand = strand;
}

void Frame::setBegin( UIntType begin )
{
    _begin = begin;
}

void Frame::setEnd( UIntType end )
{
    _end = end;
}

IdType Frame::getContigId() const
{
    return _ctgId;
}

char Frame::getStrand() const
{
    return _strand;
}

UIntType Frame::getBegin() const
{
    return _begin;
}

UIntType Frame::getEnd() const
{
    return _end;
}

UIntType Frame::getLength() const
{
    return _end - _begin + 1;
}

const Frame& 
Frame::operator=(const Frame& frame) 
{
    this->_ctgId = frame._ctgId;
    this->_begin = frame._begin;
    this->_end = frame._end;
    this->_strand = frame._strand;
    
    return *this;
}


bool
Frame::operator <(const Frame& frame) const
{
    return( this->_ctgId < frame.getContigId() || 
            (this->_ctgId == frame.getContigId() && this->_begin < frame.getBegin()) ||
            (this->_ctgId == frame.getContigId() && this->_begin == frame.getBegin() && this->_end < frame.getEnd())
            );
}


bool
Frame::operator ==(const Frame& frame) const
{
    return( this->_ctgId == frame.getContigId() &&
            this->_strand == frame.getStrand() &&
            this->_begin == frame.getBegin() &&
            this->_end == frame.getEnd() );
}


std::ostream &operator<<( std::ostream &output, const Frame &frame )
{
    output << frame.getContigId() << "\t" << frame.getStrand() << "\t" 
           << frame.getBegin() << "\t" << frame.getEnd();
    
    return output;
}

std::istream &operator>>( std::istream &input, Frame &frame )
{
    IdType ctgId;
    char strand;
    UIntType begin, end;
    
    input >> ctgId >> strand >> begin >> end;
    
    frame = Frame( ctgId, strand, begin, end);
    
    return input;
}