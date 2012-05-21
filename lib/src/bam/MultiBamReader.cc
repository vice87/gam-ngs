#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "bam/MultiBamReader.hpp"
#include "UtilityFunctions.hpp"

MultiBamReader::MultiBamReader() :
	_bam_readers(),
	_bam_aligns(),
	_valid_aligns(),
	_isize_mean(),
	_isize_std(),
	_isize_count()
{}


MultiBamReader::~MultiBamReader()
{
	this->Close();
	for( size_t i=0; i < _bam_readers.size(); i++ ) delete _bam_readers[i];
}

bool MultiBamReader::Open( const std::vector< std::string > &filenames )
{
	size_t bams = filenames.size();

	_bam_readers.resize( bams );
	_bam_aligns.resize( bams );
	_valid_aligns.resize( bams );

	_isize_mean.resize( bams, 0 );
	_isize_std.resize( bams, 0 );
	_isize_count.resize( bams, 1 );

	std::string index_filename;
	bool ret = true;

	for( size_t i=0; i < bams; i++ )
	{
		_bam_readers[i] = new BamReader();

		if( not _bam_readers[i]->Open( filenames[i] ) )
		{
			ret = false;
			std::cerr << "Unable to open BAM file:\n" << filenames[i] << std::endl;
		}
		else
		{
			std::cerr << filenames[i] << " successfully opened!" << std::endl;
		}
	}

	this->OpenIndices( filenames ); // try opening indices files (same filenames with .bai extension)

	for( size_t i=0; i < bams; i++ )
	{
		_valid_aligns[i] = _bam_readers[i]->GetNextAlignment( _bam_aligns[i] );
	}

	return ret;
}


bool MultiBamReader::Open( const std::string &filename )
{
	std::vector< std::string > filenames(1,filename);
	return this->Open( filenames );
}


bool MultiBamReader::OpenIndices( const std::vector< std::string > &filenames )
{
	std::string index_filename;
	bool ret = true;

	for( size_t i=0; i < filenames.size(); i++ )
	{
		index_filename = filenames[i] + ".bai";

		if( not _bam_readers[i]->OpenIndex( index_filename ) )
		{
			ret = false;
			std::cerr << "Unable to open BAM index file:\n" << index_filename << std::endl;
		}
	}

	return ret;
}


void MultiBamReader::Close()
{
	for( size_t i=0; i < _bam_readers.size(); i++ ) _bam_readers[i]->Close();
}


uint32_t MultiBamReader::size() const
{
	return _bam_readers.size();
}


BamReader* MultiBamReader::getBamReader( uint32_t idx )
{
	if( idx >= _bam_readers.size() ) throw MultiBamReaderException( "getBamReader index out of bound." );

	return _bam_readers[idx];
}


double MultiBamReader::getISizeMean( uint32_t idx )
{
	if( idx >= _isize_mean.size() ) throw MultiBamReaderException( "getISizeMean index out of bound." );

	return _isize_mean[idx];
}


double MultiBamReader::getISizeStd( uint32_t idx )
{
	if( idx >= _isize_std.size() ) throw MultiBamReaderException( "getISizeStd index out of bound." );

	return _isize_std[idx];
}



uint64_t MultiBamReader::getISizeNum( uint32_t idx )
{
	if( idx >= _isize_count.size() ) throw MultiBamReaderException( "getISizeNum index out of bound." );

	return _isize_count[idx];
}


bool MultiBamReader::Rewind()
{
	bool ret = true;

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		if( not _bam_readers[i]->Rewind() ) ret = false; else _valid_aligns[i] = _bam_readers[i]->GetNextAlignment( _bam_aligns[i] );
	}

	return ret;
}


bool MultiBamReader::Jump( uint32_t refID, uint32_t position )
{
	bool ret = true;

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		if( not _bam_readers[i]->Jump( refID, position ) )
		{
			ret = false;
			_valid_aligns[i] = false;
		}
		else
		{
			_valid_aligns[i] = _bam_readers[i]->GetNextAlignment( _bam_aligns[i] );
		}
	}

	return ret;
}


bool MultiBamReader::SetRegion ( const uint32_t &leftRefID, const uint32_t &leftPosition, const uint32_t &rightRefID, const uint32_t &rightPosition )
{
	bool ret = true;

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		if( not _bam_readers[i]->SetRegion( leftRefID, leftPosition, rightRefID, rightPosition ) )
		{
			ret = false;
			_valid_aligns[i] = false;
		}
		else
		{
			_valid_aligns[i] = _bam_readers[i]->GetNextAlignment( _bam_aligns[i] );
		}
	}

	return ret;
}


const RefVector& MultiBamReader::GetReferenceData() const
{
	if( _bam_readers.size() == 0 ) throw MultiBamReaderException( "GetReferenceData called on empty MultiBamReader object" );

	return (_bam_readers.front())->GetReferenceData();
}


bool MultiBamReader::GetNextAlignment( BamAlignment &align, bool update_stats )
{
	if( _bam_readers.size() == 0 ) return false;

	bool found = false;
	size_t idx = 0;

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		if( not found ) // if a valid alignment hasn't been found yet
		{
			if( not _valid_aligns[i] ) continue;

			found = true;
			align = _bam_aligns[i];
			idx = i;
		}
		else
		{
			if( _valid_aligns[i] )
			{
				if( (_bam_aligns[i].RefID == align.RefID && _bam_aligns[i].Position < align.Position) || _bam_aligns[i].RefID < align.RefID )
				{
					align = _bam_aligns[i];
					idx = i;
				}
			}
		}
	}

	// update alignments vector retrieving a new one from the proper BamReader
	if( found )
	{
		_valid_aligns[idx] = _bam_readers[idx]->GetNextAlignment( _bam_aligns[idx] );

		// update statistics
		if( update_stats && align.IsMapped() && !align.IsDuplicate() && align.IsPrimaryAlignment() && !align.IsFailedQC() &&
			align.IsFirstMate() && align.IsMateMapped() && align.RefID == align.MateRefID )
		{
			int32_t alignmentLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t startPaired = align.MatePosition;

			int32_t iSize;

			if( startRead < startPaired )
			{
				iSize = (startPaired + align.Length) - startRead;
				if( iSize < MIN_ISIZE || iSize > MAX_ISIZE ) return found;

				if( !align.IsReverseStrand() && align.IsMateReverseStrand() )
				{
					if(_isize_count[idx] == 1)
					{
						_isize_mean[idx] = iSize;
						_isize_std[idx] = 0;
						_isize_count[idx]++;
					}
					else
					{
						double oldMean = _isize_mean[idx];
						double oldStd = _isize_std[idx];

						_isize_mean[idx] = oldMean + (iSize - oldMean)/double(_isize_count[idx]);
						_isize_std[idx] = oldStd + (_isize_count[idx]-1)*(iSize - oldMean)*(iSize - oldMean)/double(_isize_count[idx]);
						_isize_count[idx]++;
					}
				}
			}
			else
			{
				iSize = (startRead + alignmentLength) - startPaired;
				if( iSize < MIN_ISIZE || iSize > MAX_ISIZE ) return found;

				if( align.IsReverseStrand() && !align.IsMateReverseStrand() )
				{
					if(_isize_count[idx] == 1)
					{
						_isize_mean[idx] = iSize;
						_isize_std[idx] = 0;
						_isize_count[idx]++;
					}
					else
					{
						double oldMean = _isize_mean[idx];
						double oldStd = _isize_std[idx];

						_isize_mean[idx] = oldMean + (iSize - oldMean)/double(_isize_count[idx]);
						_isize_std[idx] = oldStd + (_isize_count[idx]-1)*(iSize - oldMean)*(iSize - oldMean)/double(_isize_count[idx]);
						_isize_count[idx]++;
					}
				}
			}
		}
	}
	else
	{
		if( update_stats ){ for( size_t i=0; i < _isize_std.size(); i++ ) _isize_std[i] = sqrt( _isize_std[i] / double(_isize_count[i]) ); }
	}

	return found;
}


void MultiBamReader::writeStatsToFile( const std::string &filename ) const
{
	std::ofstream ofs( filename.c_str() );

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		ofs << getPathBaseName(_bam_readers[i]->GetFilename()) << std::endl;
		ofs << _isize_mean[i] << "\t" << _isize_std[i] << "\t" << _isize_count[i] << std::endl;
	}

	ofs.close();
}

uint32_t MultiBamReader::readStatsFromFile( const std::string &filename )
{
	std::ifstream ifs( filename.c_str() );

	std::string bamfile, data;
	uint32_t idx = 0;

	while( ifs.good() )
	{
		getline( ifs, bamfile ); // read filename
		if( bamfile == "" ) continue;

		if( idx >= _bam_readers.size() )
		{
			std::cerr << "The number of inserts statistics doesn't match the number of BAM files.\nStats file: " << filename << std::endl;
			return idx;
		}

		if( bamfile != getPathBaseName(_bam_readers[idx]->GetFilename()) )
		{
			std::cerr << "Error loading bam statistics file (corresponding BAM file not found): " << bamfile << std::endl;
			return idx;
		}

		getline( ifs, data );
		std::stringstream ss(data);
		ss >> _isize_mean[idx] >> _isize_std[idx] >> _isize_count[idx];

		idx++;
	}

	ifs.close();

	return idx;
}

BamReader& MultiBamReader::at( const size_t &index ) const
{
	return *(this->_bam_readers.at(index));
}


BamReader& MultiBamReader::operator[](const size_t &index ) const
{
	return *(this->_bam_readers[index]);
}