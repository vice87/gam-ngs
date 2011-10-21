/* 
 * File:   ContigMemPool.hpp
 * Author: Riccardo Vicedomini
 *
 * Created on 4 giugno 2011, 23.57
 */

#ifndef CONTIGMEMPOOL_HPP
#define	CONTIGMEMPOOL_HPP

#include "assembly/contig.hpp"

class ContigMemPool : public std::map<std::string,Contig>
{
    friend std::istream& operator >>(std::istream &is, ContigMemPool &cmp);
    friend std::ostream& operator <<(std::ostream &os, const ContigMemPool& cmp);
    
public:
    
    Contig get( const std::string &id ) const;
    void set( const std::string &id, const Contig &ctg );
    
    static ContigMemPool loadPool(const std::string &file);
    static ContigMemPool loadPool(const char *file);
    
    static void savePool(const std::string& poolFile, const ContigMemPool &pool);
    static void savePool(const char* poolFile, const ContigMemPool &pool);
};

#endif	/* CONTIGMEMPOOL_HPP */

