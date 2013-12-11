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

/*!
 * \file StrandProbability.hpp
 * \brief Definition of StrandProbability class.
 * \details This file contains the definition of the class representing a
 * probability value.
 */

#ifndef STRANDPROBABILITY_HPP
#define	STRANDPROBABILITY_HPP

#include "types.hpp"

//! Class StrandProbability.
class StrandProbability
{

private:
    RealType _probability; //!< 0 means "-", 1 means "+"

    //! Limits the probability in the range [0,1]
    void boundValue();

public:
    //! A constructor.
    /*!
     * Initialises the probability to 0.5.
     */
    StrandProbability();

    //! A copy constructor.
    /*!
     * \param orig a StrandProbability object.
     */
    StrandProbability( const StrandProbability &orig );

    //! A constructor.
    /*!
     * Initialises the object with a given probability.
     * \param probability a real number in [0,1] interval.
     */
    StrandProbability( const RealType &probability );

    //! Assign operator, given a StrandProbability object
    const StrandProbability& operator=( const StrandProbability &orig );

    //! Assign operator, given a probability.
    const StrandProbability& operator=( const RealType &probability );

    //! Product operator, given a StrandProbability object
    StrandProbability operator*( const StrandProbability &sp ) const;

    //! Product operator, given a probability.
    StrandProbability operator*( const RealType &val ) const;

    //! Returns a character representing the probability value.
    /*!
     * \return \c '+' if the probability is > 0.5, \c '-' if the probability is > 0.5, '?' otherwise.
     */
    char getStrand() const;

    //! Cast to RealType operator.
    /*!
     * \sa RealType
     */
    operator RealType() const;
};


#endif	/* STRANDPROBABILITY_HPP */

