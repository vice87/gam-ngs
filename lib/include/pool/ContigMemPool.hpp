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

