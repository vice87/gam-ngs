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

#ifndef HASHCONTIGMEMPOOL_HPP
#define	HASHCONTIGMEMPOOL_HPP

#include <map>
#include <set>

#include "api/BamAux.h"

#include "assembly/contig.hpp"

class HashContigMemPool
{
private:
    typedef std::map< std::string, Contig > ContigMap;
    typedef std::map< std::string, int32_t > RefMap;
    ContigMap _pool;

public:
    HashContigMemPool();

	uint64_t size() const;

    const Contig& get( const std::string &id ) const;
    void set( const std::string &id, const Contig &ctg );

	void getNames(std::set<std::string> &ctgNames) const;

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

class ExtContigMemPool
{
private:
    typedef std::map< std::string, int32_t > RefMap;
    std::vector< HashContigMemPool > _poolVect;

public:
    ExtContigMemPool();
    ExtContigMemPool( size_t num );

	uint64_t size() const;
	uint64_t size( size_t aId ) const;

    const Contig& get( const size_t aId, const std::string &name ) const;
    void set( const size_t aId, const std::string &name, const Contig &ctg );

	// carica i nomi dei contig nel pool (con assembly id specificato) in una struttura di tipo set
	void getNames(const size_t aId, std::set<std::string> &ctgNames) const;

    void loadPool(const size_t aId, const std::string &file, RefMap &refMap);

    void resize( size_t num );
    void clear();
};

#endif	/* HASHCONTIGMEMPOOL_HPP */

