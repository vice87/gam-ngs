#ifndef MULTI_BAM_READER_H_
#define MULTI_BAM_READER_H_

#include <string>
#include <vector>
#include <exception>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#define MIN_ISIZE 100
#define MAX_ISIZE 1000000

using namespace BamTools;

class MultiBamReaderException : public std::exception
{
	std::string _what;

public:
	MultiBamReaderException( const std::string &what ) : _what(what) {}
	MultiBamReaderException( const char* what ) : _what( std::string(what) ) {}

	~MultiBamReaderException() throw() {}

	virtual const char* what() const throw(){ return _what.c_str(); }
};


//! class that can handle multiple bam files of different libraries aligned on the same assembly
class MultiBamReader
{
private:
	std::vector< BamReader* > _bam_readers; 	// pointers to BAM readers
	std::vector< BamAlignment > _bam_aligns; 	// Next alignment to be processed for each reader
	std::vector< bool > _valid_aligns;			// Whether an alignment is valid (to be processed)

	std::vector< double > _isize_mean;			// mean insert size for the libraries
	std::vector< double > _isize_std;			// standard deviation of insert sizes for the libraries
	std::vector< uint64_t > _isize_count;

public:
	MultiBamReader();
	~MultiBamReader();

	bool Open( const std::vector< std::string > &filenames );
	bool Open( const std::string &filename );
	void Close();

	uint32_t size() const;

	BamReader* getBamReader( uint32_t idx );
	double getISizeMean( uint32_t idx );
	double getISizeStd( uint32_t idx );
	uint64_t getISizeNum( uint32_t idx );

	bool Rewind();
	bool Jump( uint32_t refID, uint32_t position = 0 );
	bool SetRegion ( const uint32_t &leftRefID, const uint32_t &leftPosition, const uint32_t &rightRefID, const uint32_t &rightPosition );

	bool GetNextAlignment( BamAlignment &align, bool update_stats = false );
	const RefVector& GetReferenceData() const;

	void writeStatsToFile( const std::string &filename ) const;
	uint32_t readStatsFromFile( const std::string &filename );

	BamReader& at( const size_t &index ) const;
	BamReader& operator[]( const size_t &index ) const;

private:
	bool OpenIndices( const std::vector< std::string > &filenames );
};

#endif /* MULTI_BAM_READER_H_ */