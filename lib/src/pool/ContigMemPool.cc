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

