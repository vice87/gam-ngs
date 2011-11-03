/* 
 * File:   ContigMemPool.hpp
 * Author: Riccardo Vicedomini
 *
 * Created on 4 giugno 2011, 23.57
 */

#ifndef HASHCONTIGMEMPOOL_HPP
#define	HASHCONTIGMEMPOOL_HPP

#include <map>
#include "assembly/contig.hpp"

class HashContigMemPool
{   
private:
    typedef std::map< std::string, Contig > ContigMap;
    typedef std::map< std::string, int32_t > RefMap;
    ContigMap _pool;
    
public:
    HashContigMemPool();
    
    Contig get( const std::string &id ) const;
    void set( const std::string &id, const Contig &ctg );
    
    void loadPool(const std::string &file, RefMap &refMap);
    void loadPool(const char *file, RefMap &refMap);
    
    void savePool(const std::string& poolFile);
    void savePool(const char* poolFile);
    
    void clear();
    
    void readNextContigID( std::istream &is, std::string &ctg_id );
    void readNextSequence( std::istream &is, Contig &ctg );
    
    //friend std::istream& operator >>(std::istream &is, HashContigMemPool &cmp);
    //friend std::ostream& operator <<(std::ostream &os, const HashContigMemPool& cmp);
};

#endif	/* HASHCONTIGMEMPOOL_HPP */

