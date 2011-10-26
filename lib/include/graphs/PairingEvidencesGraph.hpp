/*! 
 * \file PairingEvidencesGraph.hpp
 * \brief Definition of PairingEvidencesGraph class.
 * \details This class construct a PairedContigGraph whose edges are weighted
 * according to the number of blocks on two contigs.
 */

#ifndef PAIRINGEVIDENCESGRAPH_HPP
#define	PAIRINGEVIDENCESGRAPH_HPP

#include "graphs/PairedGraph.code.hpp"

//! PairingEvidencesGraph class.
/*!
 * Extends the PairedContigGraph class, adding a weight to each edge.
 */
class PairingEvidencesGraph : public PairedContigGraph<>
{
    
public:
    typedef PairedContigGraph<>::Edge Edge; //!< edge descriptor type
    typedef PairedContigGraph<>::Vertex Vertex; //!< vertex descriptor type
    
private: 
    //! Adds edge weights.
    /*!
     * For each block, the weight of the edge connecting its master and slave contigs
     * is increased by 1.
     * \param blocks a vector of blocks.
     */
    void addEdgeWeights( const std::vector<Block>& blocks );

public:
    //! A constructor.
    /*!
     * \param blocks a vector of blocks.
     */
    PairingEvidencesGraph( const std::vector<Block> &blocks );

};

//! Discards blocks whose master and slave contigs have a pairing evidence below a certain threshold.
/*!
 * \param blocks a vector of blocks.
 * \param minPairEvid minimum number of blocks between two contigs.
 * \return a filtered vector of blocks.
 */
std::vector<Block> filterBlocksByPairingEvidences( const std::vector<Block> &blocks, const int minPairEvid = 1 );

#endif	/* PAIRINGEVIDENCESGRAPH_HPP */

