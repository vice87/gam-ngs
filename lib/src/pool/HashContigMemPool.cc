
#include <map>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>

#include "pool/HashContigMemPool.hpp"
#include "assembly/io_contig.hpp"
#include "types.hpp"

#define BUFFER_LEN 16384

HashContigMemPool::HashContigMemPool()
{
    _pool.clear(); //_pool.resize(0);
}


uint64_t HashContigMemPool::size() const
{
	return (this->_pool).size();
}


const Contig& HashContigMemPool::get(const std::string& name) const
{
    ContigMap::const_iterator pos = (this->_pool).find(name);

    if( pos == (this->_pool).end() ) throw std::domain_error("The contig associated to " + name + " is not in the pool");
    return pos->second;
}

void HashContigMemPool::set(const std::string &name, const Contig &ctg)
{
    this->_pool[ name ] = ctg;
}


void HashContigMemPool::getNames(std::set<std::string> &ctgNames) const
{
	ContigMap::const_iterator seq = (this->_pool).begin();
	while( seq != (this->_pool).end() ){ ctgNames.insert(seq->first); ++seq; }
}


void HashContigMemPool::loadPool(const std::string &file, RefMap &refMap)
{
    std::ifstream ifs( file.c_str(), std::ifstream::in );

    char buffer[BUFFER_LEN];
    ifs.rdbuf()->pubsetbuf( buffer, BUFFER_LEN );

    while( !ifs.eof() )
    {


        std::string ctg_name;
        readNextContigID( ifs, ctg_name );

        Contig *ctg = &(this->_pool[ ctg_name ]);

        ctg->set_name( ctg_name );

		RefMap::const_iterator ref = refMap.find(ctg_name);
		if( ref != refMap.end() ) ctg->resize( ref->second );

        readNextSequence( ifs, *ctg );
    }

    ifs.close();
}

void HashContigMemPool::loadPool(const char *file, RefMap &refMap)
{
    return HashContigMemPool::loadPool( std::string(file), refMap );
}

void HashContigMemPool::readNextContigID( std::istream &is, std::string &ctg_id )
{
    char c = is.peek();

    while( is.good() and (c == ' ' or c == '\n') )
    {
        is.ignore(1);
        if( is.good() ) c = is.peek();
    }

    if( is.good() and c != '>' )
    {
        std::stringstream ss;
        ss << "Found invalid character: " << c;
	throw std::domain_error(ss.str().c_str());
    }

    std::string line, id;

    // get name
    std::getline( is, line );
    id = line.substr(1,line.size()-1);

    size_t pos = id.find(' ');
    if( pos != std::string::npos ) id = id.substr(0,pos);

    ctg_id = id;
}

void HashContigMemPool::readNextSequence( std::istream &is, Contig &ctg )
{
    if( is.eof() ) return;

    char c('\n');
    UIntType idx = 0;

    // while we do not reach a new contig
    while( !is.eof() and c != '>' )
    {
        // read a new char
        is.get(c);

        if( c != '\n' and c != '>' and c != ' ' and !is.eof() )
        {
            // copy read nucleotide into sequence
            if(idx >= ctg.size()) ctg.resize(idx+1);

            ctg.at(idx) = c;
            idx++;
        }
    }

    if( !is.eof() ) is.unget();
}


void HashContigMemPool::savePool(const std::string& poolFile)
{
    std::ofstream os(poolFile.c_str(),std::ios::out);

    ContigMap::const_iterator iter;
    for( iter = (this->_pool).begin(); iter != (this->_pool).end(); ++iter )
    {
        os << iter->second << std::endl;
    }

    os.close();
}

void HashContigMemPool::savePool(const char* poolFile)
{
    HashContigMemPool::savePool(std::string(poolFile));
}

void HashContigMemPool::clear()
{
    (this->_pool).clear();
}


ExtContigMemPool::ExtContigMemPool() :
        _poolVect(1)
{}

ExtContigMemPool::ExtContigMemPool(size_t num) :
        _poolVect(num)
{}

uint64_t ExtContigMemPool::size() const
{
	uint64_t sequences = 0;

	for( size_t i = 0; i < (this->_poolVect).size(); i++ )
		sequences += (this->_poolVect).at(i).size();

	return sequences;
}

uint64_t ExtContigMemPool::size( size_t aId ) const
{
	if( aId >= (this->_poolVect).size() ) return 0;

	return (this->_poolVect).at(aId).size();
}

const Contig& ExtContigMemPool::get(const size_t aId, const std::string& name) const
{
    return _poolVect[aId].get(name);
}

void ExtContigMemPool::set(const size_t aId, const std::string& name, const Contig& ctg)
{
    _poolVect[aId].set( name, ctg );
}

void ExtContigMemPool::getNames(const size_t aId, std::set<std::string> &ctgNames) const
{
	ctgNames.clear();
	_poolVect[aId].getNames(ctgNames);
}

void ExtContigMemPool::loadPool(const size_t aId, const std::string& file, RefMap& refMap)
{
    _poolVect[aId].loadPool(file,refMap);
}

void ExtContigMemPool::resize(size_t num)
{
    _poolVect.resize(num);
}

void ExtContigMemPool::clear()
{
    std::vector< HashContigMemPool >::iterator pool;

    for( pool = _poolVect.begin(); pool != _poolVect.end(); pool++ )
    {
        pool->clear();
    }
}