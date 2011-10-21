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

