
#include <iosfwd>
#include <map>
#include <stdexcept>
#include <ios>
#include <iostream>
#include <fstream>

#include "pool/HashContigMemPool.hpp"
#include "assembly/io_contig.hpp"

HashContigMemPool::HashContigMemPool()
{
    _pool.resize(0);
}

Contig HashContigMemPool::get(const std::string& name) const
{
    ContigMap::const_iterator pos = (this->_pool).find(name);
    
    if( pos == (this->_pool).end() ) throw std::domain_error("The contig associated to " + name + " is not in the pool");
    return pos->second;
}

void HashContigMemPool::set(const std::string &name, const Contig &ctg)
{
    this->_pool[ name ] = ctg;
}


void HashContigMemPool::loadPool(const std::string &file)
{
    std::ifstream ifs( file.c_str(), std::ifstream::in );
    (this->_pool).resize(0);
    
    while( !ifs.eof() )
    {
        Contig ctg;
        
        ifs >> ctg;
        this->_pool[ ctg.name() ] = ctg;
    }
    
    ifs.close();
}


void HashContigMemPool::loadPool(const char *file)
{
    return HashContigMemPool::loadPool( std::string(file) );
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