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

#ifndef _NUCLEOTIDE_
#define _NUCLEOTIDE_

//! Enum representing the type of a base.
typedef enum
{
    A, /*!< Adenine */
    T, /*!< Thymine */
    C, /*!< Cytosine */
    G, /*!< Guanine */
    N, /*!< Any of the four types. */
    LAST_BASE
} __attribute__((packed)) BaseType;

//! Class implementing a nucleotide.
class Nucleotide
{
private:
    BaseType _base; //!< base type.

public:
    Nucleotide();
    Nucleotide(const Nucleotide& orig);
    Nucleotide(const BaseType base_val);
    Nucleotide(const char base);

    const BaseType&
    base() const;

    const Nucleotide&
    operator=(const Nucleotide& orig);

    const Nucleotide&
    operator=(const char orig);

    bool
    operator==(const Nucleotide& a);

    bool
    operator!=(const Nucleotide& a);

    operator char() const;

    friend Nucleotide complement(const Nucleotide &orig);
};

Nucleotide complement(const Nucleotide &orig);
std::ostream& operator<<(std::ostream& os, const Nucleotide& base);

#endif // _NUCLEOTIDE_
