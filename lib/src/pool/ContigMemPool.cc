
#include <iosfwd>
#include <map>
#include <stdexcept>
#include <ios>
#include <iostream>
#include <fstream>

#include "pool/ContigMemPool.hpp"
#include "assembly/io_contig.hpp"

Contig ContigMemPool::get(const std::string& name) const
{
    ContigMemPool::const_iterator pos = this->find(name);
    if( pos == this->end() ) throw std::domain_error("The contig associated to "+name+" is not in the pool");
    return pos->second;
}

void ContigMemPool::set(const std::string &name, const Contig &ctg)
{
    (this->operator [](name)) = ctg;
}


ContigMemPool ContigMemPool::loadPool(const std::string &file)
{
    ContigMemPool pool;
    std::ifstream ifs( file.c_str(), std::ifstream::in );
    ifs >> pool;
    ifs.close();
    
    return pool;
}


ContigMemPool ContigMemPool::loadPool(const char *file)
{
    return ContigMemPool::loadPool( std::string(file) );
}


void ContigMemPool::savePool(const std::string& poolFile, const ContigMemPool &pool)
{
    std::ofstream out(poolFile.c_str(),std::ios::out);
    out << pool;
    out.close();
}

void ContigMemPool::savePool(const char* poolFile, const ContigMemPool &pool)
{
    ContigMemPool::savePool(std::string(poolFile),pool);
}


std::istream& operator >>(std::istream &is, ContigMemPool &cmp)
{
    cmp.clear();
    
    while( !is.eof() )
    {
        Contig ctg;
        
        is >> ctg;
        cmp[ ctg.name() ] = ctg;
    }
    
    return is;
}


std::ostream& operator <<(std::ostream &os, const ContigMemPool &cmp)
{
    ContigMemPool::const_iterator iter;
    for( iter = cmp.begin(); iter != cmp.end(); ++iter )
    {
        os << iter->second << std::endl;
    }
    
    return os;
}
		
		