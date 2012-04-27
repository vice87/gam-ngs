/*! 
 * \file Read.hpp
 * \brief Definition of Read class.
 * \details This file contains the definition of the class representing a read.
 */

#ifndef READS_HPP
#define	READS_HPP

#include <map>
#include <string>

#include "api/BamReader.h"

#include "types.hpp"
#include "google/sparse_hash_map"

using namespace BamTools;
using google::sparse_hash_map;

//! Class implementing a read.
class Read 
{

private: 
    int32_t _contigId;   //!< contig's identifier.
    int32_t _startPos;   //!< starting position (0-based) of the read inside the contig.
    int32_t _endPos;     //!< ending position (0-based) of the read inside the contig (half-open interval).
    bool _isRev;        //!< whether the read is reverse complemented or not.

public:
    //! A constructor with no arguments.
    Read();
    
    //! A copy constructor.
    /*!
     * Creates a copy of a read.
     * \param orig a Read object.
     */
    Read(const Read &orig);
    
    //! A constructor which sets the read's attributes.
    /*! 
     * \param ctg       contig's identifier
     * \param sPos      starting position (0-based) inside the contig
     * \param ePos      end position (0-based) inside the contig
     * \param rev       \c true if the read is reverse complemented
     */
    Read(const int32_t ctg, const int32_t sPos, const int32_t ePos, const bool rev = false);
    
    //! Returns contig's identifier.
    /*!
     * \return contig's identifier.
     */
    int32_t getContigId() const;
    
    //! Returns starting position of the read in the contig
    /*!
     * \return starting position (0-based)
     */
    int32_t getStartPos() const;
    
    //! Returns ending position of the read in the contig
    /*!
     * \return ending position (0-based)
     */
    int32_t getEndPos() const;
    
    //! Returns the length of the read.
    /*!
     * \return read's length
     */
    int32_t getLength() const;
    
    //! Returns whether the read is reverse complemented or not.
    /*!
     * \return \c true if the read is reverse complemented, \c false otherwise.
     */
    bool isReverse();
    
    bool overlaps( Read &read, int minOverlap = 0 ) const;
    
    //! Loads a set of common reads from two bam files sorted by read's name.
    /*!
     * In particular, exclusively the reads from the master bam are loaded.
     * Reads with multiple alignments or unmapped are discarded.
     * 
     * \param bamMaster alignment of reads in master assembly (read-name sorted).
     * \param bamSlave alignment of reads in slave assembly (read-name sorted).
     * \param readMap map where the reads are loaded.
     */
    static void loadMasterReadsMap( BamReader &bamMaster, BamReader &bamSlave, sparse_hash_map< std::string, Read > &readMap );
    
    //! Loads a set of (mapped) reads from a bam file.
    /*!
     * Reads with multiple alignments or unmapped are discarded.
     * 
     * \param bamReader BamReader object.
     * \param readMap map where the reads are loaded (output)
     * \param coverage vector of coverages of the contigs (output)
     */
    static void loadReadsMap( 
        BamReader &bamReader, 
        sparse_hash_map< std::string, Read > &readMap_1,
        sparse_hash_map< std::string, Read > &readMap_2,
        std::vector< std::vector<uint32_t> > &coverage );
};

#endif	/* READS_HPP */

