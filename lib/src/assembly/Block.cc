/* 
 * File:   Block.cc
 * Author: Riccardo Vicedomini
 * 
 * Created on 22 maggio 2011, 17.40
 */

#include <cstdlib>
#include <iostream>
#include <list>
#include <unistd.h>

#include "assembly/Block.hpp"
#include "OrderingFunctions.hpp"

Block::Block():
        _numReads(0)
{}

Block::Block(const Block &block):
        _numReads(block._numReads), 
        _masterFrame(block._masterFrame), _slaveFrame(block._slaveFrame)
{}

void Block::setReadsNumber( IntType nr )
{
    _numReads = nr;
}

void Block::setMasterFrame( Frame master )
{
    _masterFrame = master;
}

void Block::setSlaveFrame( Frame slave )
{
    _slaveFrame = slave;
}

IntType Block::getReadsNumber() const
{
    return _numReads;
}

const Frame& Block::getMasterFrame() const
{
    return _masterFrame;
}

const Frame& Block::getSlaveFrame() const
{
    return _slaveFrame;
}

bool Block::isEmpty()
{
    return (_numReads == 0);
}

void Block::clear()
{
    _numReads = 0;
}

bool Block::overlaps( Read &mRead, Read &sRead, int minOverlap )
{
    if( this->isEmpty() ) return true;          
    
    // reads must share master and slave contigs of the block.
    if( mRead.getContigId() != _masterFrame.getContigId() || sRead.getContigId() != _slaveFrame.getContigId() ) return false;
    
    IntType mStartDiff = IntType(_masterFrame.getEnd()) - IntType(mRead.getStartPos()) + 1;
    IntType mEndDiff = IntType(mRead.getEndPos()) - IntType(_masterFrame.getBegin()) + 1;
    IntType sStartDiff = IntType(_slaveFrame.getEnd()) - IntType(sRead.getStartPos()) + 1;
    IntType sEndDiff = IntType(sRead.getEndPos()) - IntType(_slaveFrame.getBegin()) + 1;
    
    return (mStartDiff >= minOverlap && mEndDiff >= minOverlap && sStartDiff >= minOverlap && sEndDiff >= minOverlap);
}

bool Block::addReads( Read &mRead, Read &sRead, int minOverlap )
{
    if( mRead.getLength() < minOverlap || sRead.getLength() < minOverlap ) return false;
    
    if( isEmpty() )
    {
        _numReads = 1;
        
        _masterFrame = Frame( mRead.getContigId(), '?', mRead.getStartPos(), mRead.getEndPos() );
        _slaveFrame  = Frame( sRead.getContigId(), '?', sRead.getStartPos(), sRead.getEndPos() );
        
        _masterFrame.setReadsLen( mRead.getLength() );
        _slaveFrame.setReadsLen( sRead.getLength() );

        return true;
    }
    
    // se la read da aggiungere si sovrappone al blocco per almeno minOverlap basi
    if( this->overlaps( mRead, sRead, minOverlap ) )
    {    
        _numReads++;
        
        _masterFrame.increaseReadsLen( mRead.getLength() );
        _slaveFrame.increaseReadsLen( sRead.getLength() );
            
        if( mRead.getStartPos() < _masterFrame.getBegin() )
        {
            _masterFrame.setStrand('-');
            _masterFrame.setBegin( mRead.getStartPos() );
        }

        if( sRead.getStartPos() < _slaveFrame.getBegin() )
        {
            _slaveFrame.setStrand('-');
            _slaveFrame.setBegin( sRead.getStartPos() );
        }

        if( mRead.getEndPos() > _masterFrame.getEnd() )
        {
            _masterFrame.setStrand('+');
            _masterFrame.setEnd( mRead.getEndPos() );
        }

        if( sRead.getEndPos() > _slaveFrame.getEnd() )
        {
            _slaveFrame.setStrand('+');
            _slaveFrame.setEnd( sRead.getEndPos() );
        }

        return true;
    }
    
    return false;
}

std::vector< Block > Block::filterBlocksByOverlaps(std::vector<Block>& blocks)
{
    std::list< Block > sortedBlocks;
    std::list< Block >::iterator block, cur, next;
    
    for( UIntType i = 0; i < blocks.size(); i++ ) sortedBlocks.push_back( blocks[i] );
    
    MasterBlocksOrderer mbo;
    sortedBlocks.sort( mbo );
    
    block = sortedBlocks.begin();
        
    while( block != sortedBlocks.end() )
    {
        cur = block;
        next = ++block;
                
        if( next != sortedBlocks.end() )
        {
            if ( Frame::frameOverlap(cur->getMasterFrame(), next->getMasterFrame()) )
            {
                sortedBlocks.erase( next );
                block = cur;
            }
        }
        else
        {
            block = next;
        }
    }
    
    SlaveBlocksOrderer sbo;
    sortedBlocks.sort( sbo );
    
    block = sortedBlocks.begin();
    while( block != sortedBlocks.end() )
    {
        cur = block;
        next = ++block;
                
        if( next != sortedBlocks.end() )
        {
            if ( Frame::frameOverlap(cur->getSlaveFrame(), next->getSlaveFrame(), 100) )
            {
                sortedBlocks.erase( next );
                block = cur;
            }
        }
        else
        {
            block = next;
        }
    }
    
    std::vector< Block > outBlocks( sortedBlocks.size() );
    
    UIntType idx = 0;
    for( block = sortedBlocks.begin(); block != sortedBlocks.end(); block++ )
    {
        outBlocks[idx] = *block;
        idx++;
    }
    
    return outBlocks;
}

std::vector< Block > Block::findBlocks( BamReader &inBamMaster, 
                                        BamReader &inBamSlave, 
                                        const int minBlockSize,
                                        sparse_hash_map< std::string, Read > &slaveReadMap)
{
    typedef sparse_hash_map< std::string, Read > ReadMap;
    typedef std::list< std::string > ReadNameList;
    typedef std::vector< Block > BlockVector;
    typedef std::list< Block > BlockList;
    
    /*    
    std::cerr << "Processing Slave's Reads... " << std::flush;
        
    // Memorizzo le read dell'assembly slave
    while( inBamSlave.GetNextAlignment(align) )
    {
        // scarto le reads che hanno molteplicità > 1
        uint32_t nh;
        if( align.GetTag(std::string("NH"),nh) && nh > 1 ) continue;
        
        // considero solo le reads che sono mappate
        if( align.IsMapped() )
        {
            Read slaveRead(IdType(align.RefID), IdType(align.Position), IdType(align.GetEndPosition()), align.IsReverseStrand());
            std::string readName = align.Name + ( (align.IsFirstMate() || !align.IsMateMapped()) ? "1" : "2");
            slaveReadMap.insert( std::make_pair(readName,slaveRead) );
        }
    }
    
    std::cerr << "done." << std::endl;
    */
    
    BamAlignment align;
    BlockVector outblocks;
    BlockList curblocks;
    std::list< ReadNameList > curreads;
        
    // Costruzione dei blocchi considerando le reads del master (ordinate per contig e posizione in esso)
    while( inBamMaster.GetNextAlignment(align) )
    {
        // scarto le reads che hanno molteplicità > 1
        uint32_t nh;
        if( align.GetTag(std::string("NH"),nh) && nh > 1 ) continue;
        
        // scarto le reads che non sono mappate
        if( align.IsMapped() )
        {            
            std::string readName = align.Name;
            Read masterRead(IdType(align.RefID), IdType(align.Position), IdType(align.GetEndPosition()), align.IsReverseStrand());
            
            if( align.IsPaired() ) readName = readName + (align.IsFirstMate() ? "1" : "2");
            
            ReadMap::iterator ref = slaveReadMap.find( readName ); // cerca la read corrente fra quelle dello slave
            
            if( ref == slaveReadMap.end() ) continue; // se la read corrente non esiste nello slave, passa alla successiva
            
            bool readsAdded = false;
            BlockList::iterator block = curblocks.begin();
            std::list< ReadNameList >::iterator readlist = curreads.begin();
            
            while( block != curblocks.end() )
            {
                if( block->addReads( masterRead, ref->second, 0 ) )
                {
                    readsAdded = true;
                    readlist->push_back(ref->first); // aggiungo un riferimento alla read in posizione corrisponende al blocco
                    
                    break;
                }
                else // if reads not added to the block
                {
                    // if block is "out of scope" save it or delete it.
                    if( block->getMasterFrame().getEnd() + 1 < masterRead.getStartPos() || block->getMasterFrame().getContigId() != masterRead.getContigId() )
                    {
                        if( block->getReadsNumber() >= minBlockSize ) outblocks.push_back( *block );
                        
                        // free allocated reads associated to the block
                        ReadNameList::iterator ril;
                        for( ril = readlist->begin(); ril != readlist->end(); ril++ )
                        {
                            slaveReadMap.erase(*ril);
                        }
                        
                        block = curblocks.erase(block);
                        readlist = curreads.erase(readlist);
                        
                        continue;
                    }
                    
                    ++block;
                    ++readlist;
                }
            }
            
            if( !readsAdded ) // crea un nuovo blocco ed aggiungilo alla lista dei blocchi correnti
            {
                Block newBlock;
                newBlock.addReads( masterRead, ref->second, 0 );
                
                ReadNameList newReadMapIterList;
                newReadMapIterList.push_back(ref->first);
                
                curblocks.push_back(newBlock);
                curreads.push_back(newReadMapIterList);
            }
        }
    }
    
    // save or delete remaining blocks
    std::list< ReadNameList >::iterator readlist = curreads.begin();
    for( BlockList::iterator block = curblocks.begin(); block != curblocks.end(); block++ )
    {
        if( block->getReadsNumber() >= minBlockSize ) outblocks.push_back( *block );
        
        ReadNameList::iterator ril;
        for( ril = readlist->begin(); ril != readlist->end(); ril++ )
        {
            slaveReadMap.erase(*ril);
        }
        
        ++readlist;
    }

    slaveReadMap.clear(); // remove remaining reads from the map
    
    return outblocks;
}


bool Block::shareMasterContig(const Block& a, const Block& b)
{
    return a.getMasterFrame().getContigId() == b.getMasterFrame().getContigId();
}


bool Block::shareSlaveContig(const Block& a, const Block& b)
{
    return a.getSlaveFrame().getContigId() == b.getSlaveFrame().getContigId();
}


bool Block::shareContig(const Block& a, const Block& b)
{
    return (shareMasterContig(a,b) && shareSlaveContig(a,b));
}


std::vector< Block > Block::readBlocks( const std::string& blockFile, int minBlockSize )
{
    std::ifstream input( blockFile.c_str() );
    
    Block block;
    std::vector< Block > blocks;
        
    while( input >> block ){ if( block.getReadsNumber() >= minBlockSize ) blocks.push_back(block); }
    
    input.close();
    
    return blocks;
}


void Block::writeBlocks( const std::string& blockFile, std::vector< Block >& blocks  )
{
    std::ofstream output( blockFile.c_str() );
    
    for( std::vector<Block>::iterator it = blocks.begin(); it != blocks.end(); it++ )
    {
        output << (*it) << std::endl;
    }
    
    output.close();
}


bool
Block::operator <(const Block& block) const
{    
    return ( 
             this->_masterFrame < block.getMasterFrame() ||
             (this->_masterFrame == block.getMasterFrame() && this->_slaveFrame < block.getSlaveFrame()) ||
             (this->_masterFrame == block.getMasterFrame() && this->_slaveFrame == block.getSlaveFrame() && this->_numReads < block.getReadsNumber())
           );
}


bool
Block::operator ==(const Block& block) const
{
    return ( this->_masterFrame == block.getMasterFrame() &&
             this->_slaveFrame == block.getSlaveFrame() &&
             this->_numReads == block.getReadsNumber());
}


std::ostream& operator <<( std::ostream &output, const Block &block )
{
    output << block.getReadsNumber() << "\t" << block.getMasterFrame() << "\t" << block.getSlaveFrame();
    
    return output;
}

std::istream& operator >>( std::istream &input, Block &block )
{
    input >> block._numReads >> block._masterFrame >> block._slaveFrame;
    
    return input;
}