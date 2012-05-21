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
#include <boost/detail/container_fwd.hpp>

#include "assembly/Block.hpp"
#include "OrderingFunctions.hpp"
#include "UtilityFunctions.hpp"

Block::Block():
        _numReads(0)
{}

Block::Block(const Block &block):
        _numReads(block._numReads), _masterFrame(block._masterFrame), _slaveFrame(block._slaveFrame)
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

UIntType Block::getReadsLen() const
{
    return (_masterFrame.getBlockReadsLen() + _slaveFrame.getBlockReadsLen()) / 2;
}

const Frame& Block::getMasterFrame() const
{
    return _masterFrame;
}

Frame& Block::getMasterFrame()
{
	return _masterFrame;
}

const Frame& Block::getSlaveFrame() const
{
    return _slaveFrame;
}

Frame& Block::getSlaveFrame()
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

    return ( _masterFrame.overlaps(mRead,minOverlap) && _slaveFrame.overlaps(sRead,minOverlap) );
}

bool Block::addReads( Read &mRead, Read &sRead, int minOverlap )
{
    if( mRead.getLength() < minOverlap || sRead.getLength() < minOverlap ) return false;

    if( isEmpty() )
    {
        _numReads = 1;

        _masterFrame = Frame( mRead.getContigId(), '?', mRead.getStartPos(), mRead.getEndPos() );
        _slaveFrame  = Frame( sRead.getContigId(), '?', sRead.getStartPos(), sRead.getEndPos() );

        _masterFrame.setBlockReadsLen( UIntType(mRead.getLength()) );
        _slaveFrame.setBlockReadsLen( UIntType(sRead.getLength()) );

        return true;
    }

    // se la read da aggiungere si sovrappone al blocco per almeno minOverlap basi
    if( this->overlaps( mRead, sRead, minOverlap ) )
    {
        _numReads++;

        //_masterFrame.increaseReadsLen( mRead.getLength() );
        //_slaveFrame.increaseReadsLen( sRead.getLength() );
        //_readsLen = _readsLen + ( mRead.getLength() != sRead.getLength() ? (mRead.getLength()+sRead.getLength())/2 : mRead.getLength() );
        _masterFrame.increaseBlockReadsLen( UIntType(mRead.getLength()) );
        _slaveFrame.increaseBlockReadsLen( UIntType(sRead.getLength()) );

        if( mRead.getStartPos() < _masterFrame.getBegin() )
        {
            //_masterFrame.setStrand('-');
            _masterFrame.setBegin( mRead.getStartPos() );
        }

        if( sRead.getStartPos() < _slaveFrame.getBegin() )
        {
            //_slaveFrame.setStrand('-');
            _slaveFrame.setBegin( sRead.getStartPos() );
        }

        if( mRead.getEndPos() > _masterFrame.getEnd() )
        {
            //_masterFrame.setStrand('+');
            _masterFrame.setEnd( mRead.getEndPos() );
        }

        if( sRead.getEndPos() > _slaveFrame.getEnd() )
        {
            //_slaveFrame.setStrand('+');
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


std::vector< Block > Block::filterBlocksByCoverage(std::vector<Block>& blocks, double t )
{
    std::vector< Block > outBlocks;

    for( std::vector<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
        double mcRatio = ((double) (b->getMasterFrame()).getBlockReadsLen()) / ((double) (b->getMasterFrame()).getReadsLen());
        double scRatio = ((double) (b->getSlaveFrame()).getBlockReadsLen()) / ((double) (b->getSlaveFrame()).getReadsLen());

        if( std::max(mcRatio,scRatio) >= t ) outBlocks.push_back(*b);
    }

    return outBlocks;
}


std::vector<Block> Block::filterBlocksByLength(
	std::vector<Block> &blocks,
	const BamTools::RefVector &masterRefVector,
	const std::vector< BamTools::RefVector > &slaveRefVectors,
	int32_t min_length )
{
	std::vector<Block> out_blocks;
	std::ofstream output( "filtered_blocks" );

	for( std::vector<Block>::const_iterator b = blocks.begin(); b != blocks.end(); b++ )
	{
		const Frame& mf = b->getMasterFrame();
		const Frame& sf = b->getSlaveFrame();

		int32_t m_len = 0.3 * masterRefVector.at(mf.getContigId()).RefLength;
		int32_t s_len = 0.3 * slaveRefVectors.at(sf.getAssemblyId()).at(sf.getContigId()).RefLength;

		int32_t m_threshold = std::min( min_length, m_len );
		int32_t s_threshold = std::min( min_length, s_len );

		if( mf.getLength() >= m_threshold || sf.getLength() >= s_threshold )
		{
			out_blocks.push_back(*b);
		}
		else
		{
			output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << std::endl;
			output << (*b) << "\n" << std::endl;
		}
	}

	output.close();
	return out_blocks;
}


void Block::findBlocks(
        std::vector< Block > &outblocks,
        MultiBamReader &bamReader,
        const int minBlockSize,
        sparse_hash_map< std::string, Read > &readsMap_1,
        sparse_hash_map< std::string, Read > &readsMap_2,
        std::vector< std::vector< uint32_t > > &coverage,
        int32_t assemblyId )
{
    typedef sparse_hash_map< std::string, Read > ReadMap;
    typedef std::list< std::string > ReadNameList;
    typedef std::vector< Block > BlockVector;

    BamAlignment align;
	std::list< Block > cur_blocks;
	std::list< std::pair<uint64_t,uint64_t> > cur_evid;
    //std::list< ReadNameList > curreads; // reads associated to a block

    // initialize slave coverage vector
    const RefVector& refVect = bamReader.GetReferenceData();
    coverage.resize( refVect.size() );
    for( uint32_t i=0; i < refVect.size(); i++ ) coverage[i].resize( refVect[i].RefLength, 0 );

    int32_t nh, xt; // molteplicità delle read (nh->standard, xt->bwa)

    // Costruzione dei blocchi considerando le reads del master (ordinate per contig e posizione in esso)
    while( bamReader.GetNextAlignment(align,true) )
    {
        // ignoro le reads che non sono mappate o di cattiva qualità
        if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

        // carico i campi di tipo stringa
        //if( !align.BuildCharData() ) continue;

        // se la molteplicità non è stata definita, assumo che sia pari ad 1
        if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard SAM format field
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // bwa field

        // scarto read con molteplicità maggiore di 1
        if( nh != 1 || xt != 'U' ) continue;

        //std::string readName = align.Name;
        //if( align.IsPaired() ) readName = readName + (align.IsFirstMate() ? "1" : "2");
        Read slaveRead(align.RefID, align.Position, align.GetEndPosition(), align.IsReverseStrand());

        // update slave vector coverage
        uint32_t read_len = align.GetEndPosition() - align.Position;
        for( int i=0; i < read_len; i++ ) coverage[align.RefID][align.Position+i] += 1;

        ReadMap::iterator ref;
        if( align.IsFirstMate() ) // cerca la read corrente fra quelle del master
        {
            ref = readsMap_1.find(align.Name);
            if( ref == readsMap_1.end() ) continue;
        }
        else
        {
            ref = readsMap_2.find(align.Name);
            if( ref == readsMap_2.end() ) continue;
        }

        // blocks extension
        bool readsAdded = false;
		std::list< Block >::iterator block = cur_blocks.begin();
		std::list< std::pair<uint64_t,uint64_t> >::iterator evid = cur_evid.begin();
        //std::list< ReadNameList >::iterator readlist = curreads.begin();

        while( block != cur_blocks.end() ) // per ogni blocco ancora in costruzione
        {
            if( block->addReads( ref->second, slaveRead ) ) // se read aggiunta al blocco
            {
                readsAdded = true;

				// update evidences of frames to be oriented in the same strand
				if( (ref->second).isReverse() == slaveRead.isReverse() ) (evid->first)++;
				else (evid->second)++;

				// readlist->push_back(ref->first); // aggiungo un riferimento alla read in posizione corrispondente al blocco
                break;
            }
            else if( block->getSlaveFrame().getEnd() + 1 < slaveRead.getStartPos() || block->getSlaveFrame().getContigId() != slaveRead.getContigId() ) // block out of scope
            {
                Frame& mf = block->getMasterFrame();
				Frame& sf = block->getSlaveFrame();

				mf.setStrand('+');
				sf.setStrand( evid->first >= evid->second ? '+' : '-' );

				if( block->getReadsNumber() >= minBlockSize ) outblocks.push_back( *block );

                // free allocated reads associated to the block
                // ReadNameList::iterator ril;
                // for( ril = readlist->begin(); ril != readlist->end(); ril++ ) masterReadMap.erase(*ril);
                // readlist = curreads.erase(readlist);

                block = cur_blocks.erase(block);
				evid = cur_evid.erase(evid);

                continue;
            }

            ++block;
			++evid;
            //++readlist;
        }

        // se la read non può essere aggiunta ad alcun blocco in costruzione, creane uno nuovo
        if( !readsAdded )
        {
			Block new_block;
			new_block.addReads( ref->second, slaveRead );
			new_block.getSlaveFrame().setAssemblyId( assemblyId );

			std::pair<uint64_t,uint64_t> pair_evid(0,0);

			//ReadNameList newReadMapIterList;
            //newReadMapIterList.push_back(ref->first);

			cur_blocks.push_back(new_block);
			cur_evid.push_back(pair_evid);
            //curreads.push_back(newReadMapIterList);
        }
    }

    // save or delete remaining blocks
    std::list< Block >::iterator block = cur_blocks.begin();
	std::list< std::pair<uint64_t,uint64_t> >::iterator evid = cur_evid.begin();
	//std::list< ReadNameList >::iterator readlist = curreads.begin();

	while( block != cur_blocks.end() )
	{
		Frame& mf = block->getMasterFrame();
		Frame& sf = block->getSlaveFrame();

		mf.setStrand('+');
		sf.setStrand( evid->first >= evid->second ? '+' : '-' );

		if( block->getReadsNumber() >= minBlockSize ) outblocks.push_back( *block );

		// free allocated reads associated to the block
		//ReadNameList::iterator ril;
		//for( ril = readlist->begin(); ril != readlist->end(); ril++ ) masterReadMap.erase(*ril);
		//++readlist;

		++block;
		++evid;
	}

	// clear the hashtables of the paired reads memorized
    readsMap_1.clear();
    readsMap_2.clear();
}


void Block::updateCoverages(
        std::vector<Block> &blocks,
        std::vector< std::vector<uint32_t> > &masterCoverage,
        std::vector< std::vector<uint32_t> > &slaveCoverage )
{
    int32_t ctgId, begin, end;
    uint64_t masterReadsLen, slaveReadsLen;
    Frame newFrame;

    for( size_t i=0; i < blocks.size(); i++ )
    {
        masterReadsLen = 0;
        slaveReadsLen = 0;

        // update master coverage

        const Frame& masterFrame = blocks[i].getMasterFrame();
        ctgId = masterFrame.getContigId();
        begin = masterFrame.getBegin();
        end = masterFrame.getEnd();

        for( size_t pos = begin; pos <= end; pos++ ) masterReadsLen += masterCoverage[ctgId][pos];

        // update slave coverage

        const Frame& slaveFrame = blocks[i].getSlaveFrame();
        ctgId = slaveFrame.getContigId();
        begin = slaveFrame.getBegin();
        end = slaveFrame.getEnd();

        for( size_t pos = begin; pos <= end; pos++ ) slaveReadsLen += slaveCoverage[ctgId][pos];

        // update current block

        newFrame = blocks[i].getMasterFrame();
        newFrame.setReadsLen( masterReadsLen );
        blocks[i].setMasterFrame( newFrame );

        newFrame = blocks[i].getSlaveFrame();
        newFrame.setReadsLen( slaveReadsLen );
        blocks[i].setSlaveFrame( newFrame );
    }
}


bool Block::shareMasterContig(const Block& a, const Block& b)
{
    return (a.getMasterFrame().getAssemblyId() == b.getMasterFrame().getAssemblyId() && a.getMasterFrame().getContigId() == b.getMasterFrame().getContigId());
}


bool Block::shareSlaveContig(const Block& a, const Block& b)
{
    return (a.getSlaveFrame().getAssemblyId() == b.getSlaveFrame().getAssemblyId() && a.getSlaveFrame().getContigId() == b.getSlaveFrame().getContigId());
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


std::vector< Block > Block::readCoreBlocks( const std::string& blockFile, int minBlockSize )
{
    std::ifstream input( blockFile.c_str() );

    Block block;
    std::vector< Block > blocks;

    while( input.good() )
    {
        int32_t ctgId;
        char strand;
        int32_t begin, end;
        IntType readsLen;
        IntType readsNum;

        input >> readsNum;

        input >> ctgId >> strand >> begin >> end >> readsLen;
        Frame mf( ctgId, strand, begin, end );

        input >> ctgId >> strand >> begin >> end >> readsLen;
        Frame sf( ctgId, strand, begin, end );

        if( readsNum >= minBlockSize )
        {
            block.setReadsNumber(readsNum);
            block.setMasterFrame(mf);
            block.setSlaveFrame(sf);

            blocks.push_back(block);
        }
    }

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


void Block::writeCoreBlocks( const std::string& blockFile, std::vector< Block >& blocks  )
{
    std::ofstream output( blockFile.c_str() );

    for( std::vector<Block>::iterator b = blocks.begin(); b != blocks.end(); b++ )
    {
        Frame mf = b->getMasterFrame();
        Frame sf = b->getSlaveFrame();

        output << b->getReadsNumber() << "\t"
               << mf.getContigId() << "\t" << mf.getStrand() << "\t" << mf.getBegin() << "\t" << mf.getEnd() << "\t" << b->getReadsLen() << "\t"
               << sf.getContigId() << "\t" << sf.getStrand() << "\t" << sf.getBegin() << "\t" << sf.getEnd() << "\t" << b->getReadsLen() << std::endl;
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
             this->_numReads == block.getReadsNumber() );
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