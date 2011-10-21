/* 
 * File:   ContigMemPool.hpp
 * Author: Riccardo Vicedomini
 *
 * Created on 4 giugno 2011, 23.57
 */

#ifndef HASHCONTIGMEMPOOL_HPP
#define	HASHCONTIGMEMPOOL_HPP

#include <google/sparse_hash_map>
#include "assembly/contig.hpp"

using google::sparse_hash_map;

class HashContigMemPool
{   
private:
    typedef sparse_hash_map< std::string, Contig > ContigMap;
    ContigMap _pool;
    
public:
    HashContigMemPool();
    
    Contig get( const std::string &id ) const;
    void set( const std::string &id, const Contig &ctg );
    
    void loadPool(const std::string &file);
    void loadPool(const char *file);
    
    void savePool(const std::string& poolFile);
    void savePool(const char* poolFile);
    
    //friend std::istream& operator >>(std::istream &is, HashContigMemPool &cmp);
    //friend std::ostream& operator <<(std::ostream &os, const HashContigMemPool& cmp);
};

#endif	/* HASHCONTIGMEMPOOL_HPP */

