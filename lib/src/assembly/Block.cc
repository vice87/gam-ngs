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
#include <fstream>
#include <list>
#include <unistd.h>
#include <boost/detail/container_fwd.hpp>

#include "assembly/Block.hpp"
#include "OrderingFunctions.hpp"
#include "UtilityFunctions.hpp"


Block::Block():
        _numReads(0)
{}

Block::Block( Read &mRead, Read &sRead, int minOverlap ):
	_numReads(0)
{
	if( mRead.getLength() < minOverlap || sRead.getLength() < minOverlap ) return;

	_numReads = 1;

	_masterFrame = Frame( mRead.getContigId(), '?', mRead.getStartPos(), mRead.getEndPos() );
	_slaveFrame  = Frame( sRead.getContigId(), '?', sRead.getStartPos(), sRead.getEndPos() );

	if( mRead.getLength() > 0 ) _masterFrame.setBlockReadsLen( UIntType(mRead.getLength()) ); else _masterFrame.setBlockReadsLen(0);
	if( sRead.getLength() > 0 ) _slaveFrame.setBlockReadsLen( UIntType(sRead.getLength()) ); else _slaveFrame.setBlockReadsLen(0);
}

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

int32_t Block::getMasterId() const
{
	return _masterFrame.getContigId();
}

const Frame& Block::getSlaveFrame() const
{
    return _slaveFrame;
}

Frame& Block::getSlaveFrame()
{
	return _slaveFrame;
}

int32_t Block::getSlaveId() const
{
	return _slaveFrame.getContigId();
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

std::vector< Block > Block::filterBlocksByOverlaps(std::list<Block>& blocks)
{
    std::list< Block > sortedBlocks;
    std::list< Block >::iterator block, cur, next;

	for( std::list< Block >::iterator b = blocks.begin(); b != blocks.end(); b++ ) sortedBlocks.push_back(*b);

    MasterBlocksOrderer mbo;
    sortedBlocks.sort( mbo );

    block = sortedBlocks.begin();

    while( block != sortedBlocks.end() )
    {
        cur = block;
        next = ++block;

        if( next != sortedBlocks.end() )
        {
            if( Frame::frameOverlap(cur->getMasterFrame(), next->getMasterFrame(), 95.0) )
            {
                std::cerr << "[debug] master overlap >=95%: (" << cur->getMasterId() << "," << cur->getSlaveId() << ") - ("
					<< next->getMasterId() << "," << next->getSlaveId() << ")" << std::endl;
				std::cerr << *cur << std::endl;
				std::cerr << *next << std::endl;
				//sortedBlocks.erase( next );
                //block = cur;
            }
        }

        block = next;
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
            if ( Frame::frameOverlap(cur->getSlaveFrame(), next->getSlaveFrame(), 95.0) )
            {
				std::cerr << "[debug] slave overlap >=95%: (" << cur->getMasterId() << "," << cur->getSlaveId() << ") - ("
					<< next->getMasterId() << "," << next->getSlaveId() << ")" << std::endl;
				std::cerr << *cur << std::endl;
				std::cerr << *next << std::endl;
                //sortedBlocks.erase( next );
                //block = cur;
            }
        }

        block = next;
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


void Block::filterBlocksByCoverage(
	std::list<Block>& blocks,
	const std::set< std::pair<int32_t,int32_t> > &slb,
	double min_cov,
	double t )
{
	//std::cerr << "low coverage block threshold = " << min_cov << "X" << std::endl;
	std::list<Block>::iterator b = blocks.begin();

	while( b != blocks.end() )
	{
		double mcRatio = ((double) (b->getMasterFrame()).getBlockReadsLen()) / ((double) (b->getMasterFrame()).getReadsLen());
		double scRatio = ((double) (b->getSlaveFrame()).getBlockReadsLen()) / ((double) (b->getSlaveFrame()).getReadsLen());

		if( std::max(mcRatio,scRatio) < t )
		{
			b = blocks.erase(b); continue;
		}
		else // FILTRAGGIO BASATO SULLA COPERTURA DEL BLOCCO RISPETTO ALLA COPERTURA MEDIA //TODO:verificare di non togliere blocchi importanti
		{
			int32_t mid = b->getMasterId();
			int32_t sid = b->getSlaveId();

			if( slb.find( std::make_pair(mid,sid) ) == slb.end() )
			{
				double m_cov = ((double) (b->getMasterFrame()).getBlockReadsLen()) / (b->getMasterFrame()).getLength();
				double s_cov = ((double) (b->getSlaveFrame()).getBlockReadsLen()) / (b->getSlaveFrame()).getLength();
				double cov = mcRatio >= scRatio ? m_cov : s_cov;

				if( cov < min_cov )
				{
					//std::cerr << "low coverage block removed: (" << b->getMasterId() << "," << b->getSlaveId() << ")" << std::endl;
					b = blocks.erase(b); continue;
				}
			}
		}

		++b;
	}
}


void Block::filterBlocksByLength(
	std::list<Block> &blocks,
	const RefSequence &masterRef,
	const RefSequence &slaveRef,
	const std::set< std::pair<int32_t,int32_t> > &slb,
	int32_t min_length )
{
	std::ofstream output( "filtered_blocks" );
	std::list<Block>::iterator b, cur, prev, next;

	MasterBlocksOrderer mbo;
	blocks.sort( mbo );

	b = blocks.begin();
	while( b != blocks.end() )
	{
		cur = b;
		next = ++b;

		const Frame& mf = cur->getMasterFrame();
		const Frame& sf = cur->getSlaveFrame();

		int32_t mid = mf.getContigId();
		int32_t sid = sf.getContigId();

		if( slb.find( std::make_pair(mid,sid) ) != slb.end() ){ prev = cur; continue; }

		double m_len = 0.3 * masterRef[mid].RefLength;
		double s_len = 0.3 * slaveRef[sid].RefLength;

		double m_threshold = std::min( double(min_length), m_len );
		double s_threshold = std::min( double(min_length), s_len );

		if( mf.getLength() < m_threshold && sf.getLength() < s_threshold )
		{
			if( cur == blocks.begin() )
			{
				if( next != blocks.end() && Frame::frameOverlap(next->getMasterFrame(), cur->getMasterFrame(), 0.5) )
				{
					output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << "\n" << *cur << "\n" << std::endl;

					b = blocks.erase(cur);
					continue;
				}
			}
			else
			{
				if( Frame::frameOverlap(prev->getMasterFrame(), cur->getMasterFrame(), 0.5) )
				{
					output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << "\n" << *cur << "\n" << std::endl;

					b = blocks.erase(cur);
					continue;
				}

				if( next != blocks.end() && Frame::frameOverlap(next->getMasterFrame(), cur->getMasterFrame(), 0.5) )
				{
					output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << "\n" << *cur << "\n" << std::endl;

					b = blocks.erase(cur);
					continue;
				}
			}
		}

		prev = cur;
	}

	SlaveBlocksOrderer sbo;
	blocks.sort( sbo );

	b = blocks.begin();
	while( b != blocks.end() )
	{
		cur = b;
		next = ++b;

		const Frame& mf = cur->getMasterFrame();
		const Frame& sf = cur->getSlaveFrame();

		int32_t mid = mf.getContigId();
		int32_t sid = sf.getContigId();

		if( slb.find( std::make_pair(mid,sid) ) != slb.end() ){ prev = cur; continue; }

		double m_len = 0.3 * masterRef[mid].RefLength;
		double s_len = 0.3 * slaveRef[sid].RefLength;

		double m_threshold = std::min( double(min_length), m_len );
		double s_threshold = std::min( double(min_length), s_len );

		if( mf.getLength() < m_threshold && sf.getLength() < s_threshold )
		{
			if( cur == blocks.begin() )
			{
				if( next != blocks.end() && Frame::frameOverlap(next->getSlaveFrame(), cur->getSlaveFrame(), 0.5) )
				{
					output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << "\n" << *cur << "\n" << std::endl;

					b = blocks.erase(cur);
					continue;
				}
			}
			else
			{
				if( Frame::frameOverlap(prev->getSlaveFrame(), cur->getSlaveFrame(), 0.5) )
				{
					output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << "\n" << *cur << "\n" << std::endl;

					b = blocks.erase(cur);
					continue;
				}

				if( next != blocks.end() && Frame::frameOverlap(next->getSlaveFrame(), cur->getSlaveFrame(), 0.5) )
				{
					output << mf.getLength() << " < " << m_threshold << " and " << sf.getLength() << " < " << s_threshold << "\n" << *cur << "\n" << std::endl;

					b = blocks.erase(cur);
					continue;
				}
			}
		}

		prev = cur;
	}

	output.close();
}


/*
void Block::filterBlocksByLength(
	std::list<Block> &blocks,
	const RefSequence &masterRef,
	const RefSequence &slaveRef,
	const std::set< std::pair<int32_t,int32_t> > &slb,
	int32_t min_length )
{
	std::list<Block>::iterator b = blocks.begin();

	while( b != blocks.end() )
	{
		const Frame& mf = b->getMasterFrame();
		const Frame& sf = b->getSlaveFrame();

		int32_t m_len = 0.3 * masterRef[mf.getContigId()].RefLength;
		int32_t s_len = 0.3 * slaveRef[sf.getContigId()].RefLength;

		int32_t m_threshold = std::min( min_length, m_len );
		int32_t s_threshold = std::min( min_length, s_len );

		if( mf.getLength() < m_threshold && sf.getLength() < s_threshold )
		{
			b = blocks.erase(b);
			continue;
		}

		++b;
	}
}
*/


void Block::findBlocks(
        std::vector< Block > &outblocks,
        MultiBamReader &bamReader,
        const int minBlockSize,
        sparse_hash_map< std::string, Read > &readsMap_1,
        sparse_hash_map< std::string, Read > &readsMap_2,
        std::vector< std::vector< uint32_t > > &coverage )
{
    typedef sparse_hash_map< std::string, Read > ReadMap;
    typedef std::list< std::string > ReadNameList;
    typedef std::vector< Block > BlockVector;

    BamAlignment align;
	std::list< Block > cur_blocks;
	std::list< std::pair<uint64_t,uint64_t> > cur_evid;

    // initialize slave coverage vector
    const RefVector& refVect = bamReader.GetReferenceData();
    coverage.resize( refVect.size() );
    for( uint32_t i=0; i < refVect.size(); i++ ) coverage[i].resize( refVect[i].RefLength, 0 );

    int32_t nh, xt; // molteplicità delle read (nh->standard, xt->bwa)

    // process reads to build blocks (updating inserts statistics) by coordinate order
    while( bamReader.GetNextAlignment(align,true) )
    {
		// skip unmapped or bad-quality reads
		if( !align.IsMapped() || align.Position < 0 || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

        // load read's moltiplicity (if the field is missing, assume it as uniquely mapped)
        if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard SAM format field
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // bwa field
		bool uniqMapRead = (nh == 1 && xt == 'U');

		if( !uniqMapRead ) continue; // skip reads mapped in multiple positions

        Read slaveRead(align.RefID, align.Position, align.GetEndPosition(), align.IsReverseStrand());

        // update slave vector coverage
        uint32_t read_len = align.GetEndPosition() - align.Position;
        for( int i=0; i < read_len; i++ ) coverage[align.RefID][align.Position+i] += 1;

		ReadMap::iterator ref;

		if( !align.IsPaired() || align.IsFirstMate() ) // find the read in the first map
		{
			ref = readsMap_1.find(align.Name);
			if( ref == readsMap_1.end() ) continue; // skip the read if it has not been mapped on the other assembly
		}
		else // if second mate, find the read in the second map
		{
			ref = readsMap_2.find(align.Name);
			if( ref == readsMap_2.end() ) continue; // skip the read if it has not been mapped on the other assembly
		}

        // try to extend one of the memorized blocks
        bool readsAdded = false;
		std::list< Block >::iterator block = cur_blocks.begin();
		std::list< std::pair<uint64_t,uint64_t> >::iterator evid = cur_evid.begin();

        while( block != cur_blocks.end() ) // for each memorized block
        {
			if( block->addReads( ref->second, slaveRead ) ) // if read has been succesfully added to the current block
            {
                readsAdded = true;

				// update evidences of frames to be oriented in the same strand
				if( (ref->second).isReverse() == slaveRead.isReverse() ) (evid->first)++;
				else (evid->second)++;

				// readlist->push_back(ref->first); // aggiungo un riferimento alla read in posizione corrispondente al blocco
                break;
            }

            bool blockOutOfScope = block->getSlaveFrame().getEnd() + 1 < slaveRead.getStartPos() ||
				block->getSlaveFrame().getContigId() < slaveRead.getContigId();

			if( !readsAdded && blockOutOfScope ) // block out of scope
            {
                Frame& mf = block->getMasterFrame();
				Frame& sf = block->getSlaveFrame();

				// set frames' strand according to the numbero of concordant/discordant reads in the block
				mf.setStrand('+');
				sf.setStrand( evid->first >= evid->second ? '+' : '-' );

				if( block->getReadsNumber() >= minBlockSize ) outblocks.push_back( *block );

                // remove block and its strand evidences
				block = cur_blocks.erase(block);
				evid = cur_evid.erase(evid);

                continue;
            }

            ++block;
			++evid;
        }

        // if the read has not been added to any existing block, create a new block.
        if( !readsAdded )
        {
			Block new_block( ref->second, slaveRead, minBlockSize );
			std::pair<uint64_t,uint64_t> strand_evid(0,0);

			cur_blocks.push_back(new_block);
			cur_evid.push_back(strand_evid);
		}
    }

    // after all reads have been processed, save or delete remaining blocks
    std::list< Block >::iterator block = cur_blocks.begin();
	std::list< std::pair<uint64_t,uint64_t> >::iterator evid = cur_evid.begin();

	while( block != cur_blocks.end() )
	{
		Frame& mf = block->getMasterFrame();
		Frame& sf = block->getSlaveFrame();

		// set frames' strand according to the numbero of concordant/discordant reads in the block
		mf.setStrand('+');
		sf.setStrand( evid->first >= evid->second ? '+' : '-' );

		if( block->getReadsNumber() >= minBlockSize ) outblocks.push_back( *block );

		// remove block and its strand evidences
		block = cur_blocks.erase(block);
		evid = cur_evid.erase(evid);
	}

    readsMap_1.clear();
    readsMap_2.clear();
}


void Block::updateCoverages(
        std::vector<Block> &blocks,
        const std::vector< std::vector<uint32_t> > &masterCoverage,
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

        for( size_t pos = begin; pos <= end; pos++ ) masterReadsLen += masterCoverage.at(ctgId).at(pos);

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
    return a.getMasterId() == b.getMasterId();
}


bool Block::shareSlaveContig(const Block& a, const Block& b)
{
    return a.getSlaveId() == b.getSlaveId();
}


bool Block::shareContig(const Block& a, const Block& b)
{
    return (shareMasterContig(a,b) && shareSlaveContig(a,b));
}


void Block::loadBlocks( const std::string& blockFile, std::list<Block> &blocks, int minBlockSize )
{
    std::ifstream ifs( blockFile.c_str() );

	Block block;
	std::string line;

	while( ifs.good() )
	{
		// discard heading line
		getline( ifs, line );

		if( line == "" ) continue; // empty line
		if( line[0] == '#' ) continue; // heading/comment line

		std::stringstream ss(line);
		if( ss >> block && block.getReadsNumber() >= minBlockSize ) blocks.push_back(block);
	}


    ifs.close();
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


void Block::writeBlocks( const std::string &blockFile, std::vector< Block > &blocks  )
{
    std::ofstream output( blockFile.c_str() );

	// block file header
	output << "# MasterAssemblyID\tMasterContigID\tMasterStrand\tMasterBegin\tMasterEnd\tMasterBlockReadsLength\tMasterReadsLength\t"
		   << "SlaveAssemblyID\tSlaveContigID\tSlaveStrand\tSlaveBegin\tSlaveEnd\tSlaveBlockReadsLength\tSlaveReadsLength\n";

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


void getNoBlocksContigs(
	const RefSequence &masterRef,
	const RefSequence &slaveRef,
	const std::list<Block> &blocks,
	std::set< int32_t > &masterNBC,
	std::set< int32_t > &slaveNBC )
{
	std::vector<bool> master_ctgWithBlock( masterRef.size(), false );
	std::vector<bool> slave_ctgWithBlock( slaveRef.size(), false );

	int32_t mid, sid;
	int32_t m_size = static_cast<int32_t>(masterRef.size());
	int32_t s_size = static_cast<int32_t>(slaveRef.size());

	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); ++b )
	{
		mid = b->getMasterId();
		sid = b->getSlaveId();

		if( mid < m_size && mid >= 0 )
		{
			master_ctgWithBlock[mid] = true;
		}
		else
		{
			std::cerr << "[getNoBlocksContigs] error: found a block with master id "
					<< mid << " when the admissible range is [0," << m_size << ")." << std::endl;

			exit(EXIT_FAILURE);
		}

		if( sid < s_size && sid >= 0 )
		{
			slave_ctgWithBlock[mid] = true;
		}
		else
		{
			std::cerr << "[getNoBlocksContigs] error: found a block with slave id "
					<< sid << " when the admissible range is [0," << s_size << ")." << std::endl;

			exit(EXIT_FAILURE);
		}
	}

	for( int32_t i=0; i < master_ctgWithBlock.size(); i++ )
		if( not master_ctgWithBlock[i] ) masterNBC.insert(i);

	for( int32_t i=0; i < slave_ctgWithBlock.size(); i++ )
		if( not slave_ctgWithBlock[i] ) slaveNBC.insert(i);
}


void getNoBlocksAfterFilterContigs(
	const RefSequence &masterRef,
	const RefSequence &slaveRef,
	const std::list<Block> &blocks,
	const std::set< int32_t > &masterNBC,
	const std::set< int32_t > &slaveNBC,
	std::set< int32_t > &masterNBC_AF,
	std::set< int32_t > &slaveNBC_AF )
{
	std::vector<bool> master_ctgWithBlockAF( masterRef.size(), false );
	std::vector<bool> slave_ctgWithBlockAF( slaveRef.size(), false );

	for( std::list<Block>::const_iterator b = blocks.begin(); b != blocks.end(); ++b )
	{
		master_ctgWithBlockAF.at( b->getMasterId() ) = true;
		slave_ctgWithBlockAF.at( b->getSlaveId() ) = true;
	}

	for( std::set< int32_t >::const_iterator it = masterNBC.begin(); it != masterNBC.end(); it++ )
		master_ctgWithBlockAF.at(*it) = true;

	for( std::set< int32_t >::const_iterator it = slaveNBC.begin(); it != slaveNBC.end(); it++ )
		slave_ctgWithBlockAF.at(*it) = true;

	for( int32_t i=0; i < master_ctgWithBlockAF.size(); i++ )
		if( not master_ctgWithBlockAF.at(i) ) masterNBC_AF.insert(i);

	for( int32_t i=0; i < slave_ctgWithBlockAF.size(); i++ )
		if( not slave_ctgWithBlockAF.at(i) ) slaveNBC_AF.insert(i);
}
