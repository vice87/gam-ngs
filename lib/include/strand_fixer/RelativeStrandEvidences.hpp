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
 * \file RelativeStrandEvidences.hpp
 * \brief Definition of RelativeStrandEvidences class.
 * \details This file contains the definition of the class representing the
 * probability of a contig/sequence to be reverse complemented or not.
 */

#ifndef RELATIVESTRANDEVIDENCES_HPP
#define	RELATIVESTRANDEVIDENCES_HPP

#include <map>
#include <iostream>
#include <exception>

#include "types.hpp"

//! Class implementing the probability of a contig to be reverse complemented or not.
class RelativeStrandEvidences
{

private:
    UIntType _positive; //!< Positive evidences.
    UIntType _negative; //!< Negative evidences.

public:
    //! A constructor.
    /*!
     * Creates an empty an object with no evidences.
     */
    RelativeStrandEvidences();

    //! A copy constructor.
    /*!
     * \param orig a RelativeStrandEvidences object
     */
    RelativeStrandEvidences(const RelativeStrandEvidences &orig);

    //! Gets positive evidences number.
    /*!
     * \return the number of positive evidences.
     */
    const UIntType& getPositiveEvidences() const;

    //! Gets negative evidences number.
    /*!
     * \return the number of negative evidences.
     */
    const UIntType& getNegativeEvidences() const;

    //! Gets the number of all evidences.
    /*!
     * \return the number of all evidences.
     */
    UIntType getEvidences() const;

    //! Increase the number of positive evidences.
    /*!
     * \param positive number of positive evidences to be added.
     */
    void addPositiveEvidences(const UIntType& positive);

    //! Increase the number of negative evidences.
    /*!
     * \param positive number of negative evidences to be added.
     */
    void addNegativeEvidences(const UIntType& negative);

    //! Gets a char rapresenting the probability of a contig to be reverse complemented or not.
    /*!
     * \return '+' if there are more positive evidences, '-' if there are more negative evidences, '?' otherwise.
     */
    char getStrand() const;

    //! Gets the probability of not being reverse complemented.
    inline RealType getPositiveStrandProb() const
    {
        return ((RealType)this->_positive) / ((RealType)this->getEvidences());
    }

    //! Gets the probability of being reverse complemented.
    inline RealType getNegativeStrandProb() const
    {
        return ((RealType)this->_negative) / ((RealType)this->getEvidences());
    }

    //! Gets the probability
    /*!
     * \param strand \c '+', \c '-' or \c '?'.
     */
    inline RealType probabilityOf( const char &strand ) const
    {
        switch(strand)
        {
            case '+': return getPositiveStrandProb();
            case '-': return getNegativeStrandProb();
            case '?': return ((RealType)1) / ((RealType)2);

            default: throw std::exception();
        }
    }

    //! Assign operator.
    const RelativeStrandEvidences& operator =(const RelativeStrandEvidences &orig);

    //! Sum operator.
    RelativeStrandEvidences operator +(const RelativeStrandEvidences &orig);

    //! Sum-assign operator.
    const RelativeStrandEvidences& operator +=(const RelativeStrandEvidences &orig);
};

#endif	/* RELATIVESTRANDEVIDENCES_HPP */

