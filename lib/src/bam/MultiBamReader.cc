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
	_isize_count(),
	_asm_size(0),
	_reads_len(),
	_coverage()
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

	_bam_mutex.resize( bams );

	_minInsert.resize( bams );
	_maxInsert.resize( bams );

	_isize_mean.resize( bams, 0 );
	_isize_std.resize( bams, 0 );
	_isize_count.resize( bams, 1 );
	
	_asm_size = 0;
	_reads_len.resize( bams, 0 );
	_coverage.resize( bams, 0 );

	std::string index_filename;
	bool ret = true;

	for( size_t i=0; i < bams; i++ )
	{
		_bam_readers[i] = new BamReader();

		if( not _bam_readers[i]->Open( filenames[i] ) )
		{
			ret = false;
			std::cerr << "[bam] warning: unable to open BAM file:\n" << filenames[i] << std::endl;
		}
		else
		{
			std::cerr << "[bam] " << filenames[i] << " successfully opened!" << std::endl;
		}

		pthread_mutex_init( &(this->_bam_mutex[i]), NULL );
	}

	this->OpenIndices( filenames ); // try opening indices files (same filenames with .bai extension)

	// initialization of min/max inserts sizes
	for( size_t i=0; i < bams; i++ ) _minInsert[i] = MIN_ISIZE;
	for( size_t i=0; i < bams; i++ ) _maxInsert[i] = MAX_ISIZE;

	// get first alignment from each bam files
	for( size_t i=0; i < bams; i++ ) _valid_aligns[i] = _bam_readers[i]->GetNextAlignment( _bam_aligns[i] );
	
	this->computeAssemblySize();

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
			std::cerr << "[bam] warning: unable to open BAM index file:\n" << index_filename << std::endl;
		}
	}

	return ret;
}


void MultiBamReader::Close()
{
	for( size_t i=0; i < _bam_readers.size(); i++ ) _bam_readers[i]->Close();
}


void MultiBamReader::setMinMaxInsertSizes( const std::vector<int32_t> &minInsert, const std::vector<int32_t> &maxInsert )
{
	if( minInsert.size() != maxInsert.size() || minInsert.size() != _bam_readers.size() )
	{
		std::cerr << "min/max insert size has not been set, previous assigned values will be used (or default)." << std::endl;
		return;
	}

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		if( _minInsert[i] > 0 && _maxInsert[i] > 0 )
		{
			_minInsert[i] = minInsert[i];
			_maxInsert[i] = maxInsert[i];
		}
		else
		{
			_minInsert[i] = MIN_ISIZE;
			_maxInsert[i] = MAX_ISIZE;
		}
	}
}


uint32_t MultiBamReader::size() const
{
	return _bam_readers.size();
}


BamReader* MultiBamReader::getBamReader( uint32_t idx )
{
	if( idx >= _bam_readers.size() ) throw MultiBamReaderException( "MultiBamReader::getBamReader index out of bound." );

	return _bam_readers[idx];
}


double MultiBamReader::getISizeMean( uint32_t idx )
{
	if( idx >= _isize_mean.size() ) throw MultiBamReaderException( "MultiBamReader::getISizeMean index out of bound." );

	return _isize_mean[idx];
}


double MultiBamReader::getISizeStd( uint32_t idx )
{
	if( idx >= _isize_std.size() ) throw MultiBamReaderException( "MultiBamReader::getISizeStd index out of bound." );

	return _isize_std[idx];
}


uint64_t MultiBamReader::getISizeNum( uint32_t idx )
{
	if( idx >= _isize_count.size() ) throw MultiBamReaderException( "MultiBamReader::getISizeNum index out of bound." );

	return _isize_count[idx];
}


double MultiBamReader::getCoverage( uint32_t idx )
{
	if( idx >= _coverage.size() ) throw MultiBamReaderException( "MultiBamReader::getISizeStd index out of bound." );
	
	return _coverage[idx];
}


double MultiBamReader::getMeanCoverage()
{
    double mean_coverage = 0;
    
    for( size_t i=0; i < _coverage.size(); i++ ) mean_coverage += _coverage[i];
    
	return (mean_coverage / _coverage.size());
}


double MultiBamReader::getGlobCoverage()
{
    double glob_coverage = 0;
    for( size_t i=0; i < _coverage.size(); i++ ) glob_coverage += _coverage[i];
    
	return glob_coverage;
}


int MultiBamReader::getProperLibrary( uint64_t gap )
{
	int idx = -1;

	for( size_t i=0; i < _isize_mean.size(); i++ )
	{
		if( _isize_mean[i] != 0 && _isize_mean[i] + _isize_std[i] >= gap ) return i;
	}

	return idx;
}


void MultiBamReader::lockBamReader( uint32_t idx )
{
	pthread_mutex_lock(&(this->_bam_mutex[idx]));
}

void MultiBamReader::unlockBamReader( uint32_t idx )
{
	pthread_mutex_unlock(&(this->_bam_mutex[idx]));
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
	if( _bam_readers.size() == 0 ) throw MultiBamReaderException( "MultiBamReader::GetReferenceData called on empty MultiBamReader object" );

	return (_bam_readers.front())->GetReferenceData();
}


bool MultiBamReader::GetNextAlignment( BamAlignment &align, bool update_stats )
{
	if( _bam_readers.size() == 0 ) return false;

	bool found = false;
	size_t libId = 0;

	// retrieve next read
	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		if( not found ) // if a valid alignment hasn't been found yet
		{
			if( not _valid_aligns[i] ) continue;

			found = true;
			align = _bam_aligns[i];
			libId = i;
		}
		else
		{
			if( _valid_aligns[i] )
			{
				if( (_bam_aligns[i].RefID == align.RefID && _bam_aligns[i].Position < align.Position) || _bam_aligns[i].RefID < align.RefID )
				{
					align = _bam_aligns[i];
					libId = i;
				}
			}
		}
	}

	// update alignments vector retrieving a new one from the proper BamReader
	if( found )
	{
		// load the read following the one extracted
		_valid_aligns[libId] = _bam_readers[libId]->GetNextAlignment( _bam_aligns[libId] );
		
		if( update_stats && align.IsMapped() && !align.IsDuplicate() && align.IsPrimaryAlignment() && !align.IsFailedQC() )
		{
			_reads_len[libId] += (align.GetEndPosition() - align.Position);
		}

		// if needed, update statistics only if the read extracted has good quality and its mate is mapped on the same contig
		if( update_stats && align.IsMapped() && !align.IsDuplicate() && align.IsPrimaryAlignment() && !align.IsFailedQC() &&
			align.IsFirstMate() && align.IsMateMapped() && align.RefID == align.MateRefID )
		{
			int32_t alignmentLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t startMate = align.MatePosition;

			int32_t iSize;

			if( startRead < startMate )
			{
				iSize = (startMate + align.Length) - startRead;
				if( iSize < _minInsert[libId] || iSize > _maxInsert[libId] ) return found;

				// if the read and its mate are properly oriented update mean and std
				if( !align.IsReverseStrand() && align.IsMateReverseStrand() )
				{
					if(_isize_count[libId] == 1)
					{
						_isize_mean[libId] = iSize;
						_isize_std[libId] = 0;
						_isize_count[libId]++;
					}
					else
					{
						double oldMean = _isize_mean[libId];
						double oldStd = _isize_std[libId];

						_isize_mean[libId] = oldMean + (iSize - oldMean)/double(_isize_count[libId]);
						_isize_std[libId] = oldStd + (_isize_count[libId]-1)*(iSize - oldMean)*(iSize - oldMean)/double(_isize_count[libId]);
						_isize_count[libId]++;
					}
				}
			}
			else
			{
				iSize = (startRead + alignmentLength) - startMate;
				if( iSize < _minInsert[libId] || iSize > _maxInsert[libId] ) return found;

				// if the read and its mate are properly oriented update mean and std
				if( align.IsReverseStrand() && !align.IsMateReverseStrand() )
				{
					if(_isize_count[libId] == 1)
					{
						_isize_mean[libId] = iSize;
						_isize_std[libId] = 0;
						_isize_count[libId]++;
					}
					else
					{
						double oldMean = _isize_mean[libId];
						double oldStd = _isize_std[libId];

						_isize_mean[libId] = oldMean + (iSize - oldMean)/double(_isize_count[libId]);
						_isize_std[libId] = oldStd + (_isize_count[libId]-1)*(iSize - oldMean)*(iSize - oldMean)/double(_isize_count[libId]);
						_isize_count[libId]++;
					}
				}
			}
		}
	}
	else // the end of all bam files has been reached
	{
		// if needed, compute standard deviation
		if( update_stats )
		{ 
			for( size_t i=0; i < _isize_std.size(); i++ )
			{
				_isize_std[i] = sqrt( _isize_std[i] / double(_isize_count[i]) );
				_coverage[libId] = (_asm_size != 0) ? _reads_len[libId] / ((double)_asm_size) : 0.0;
			}
		}
	}

	return found;
}


uint64_t MultiBamReader::computeAssemblySize()
{
	_asm_size = 0;
	
	const RefVector& ref_vect = this->GetReferenceData();
	for( size_t i=0; i < ref_vect.size(); i++ ) _asm_size += ref_vect[i].RefLength;
	
	return _asm_size;
}


bool MultiBamReader::computeStatistics()
{
	if( (this->_bam_readers).size() == 0 ) return false;
	
	BamAlignment align;
	
	// for each library compute its statistics
	for( size_t libId=0; libId < _bam_readers.size(); libId++ )
	{
		_bam_readers[libId]->Rewind();
		
		this->_isize_mean[libId] = 0;
		this->_isize_std[libId] = 0;
		this->_isize_count[libId] = 1;
		
		this->_reads_len[libId] = 0;
		
		while( _bam_readers[libId]->GetNextAlignmentCore(align) )
		{
			// skip unmapped or bad-quality reads
			if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;
			
			int32_t alignmentLength = align.GetEndPosition() - align.Position;
			int32_t startRead = align.Position;
			int32_t startMate = align.MatePosition;
			
			// update reads' length
			this->_reads_len[libId] += alignmentLength;
			
			// update insert statistics only if the read extracted has its mate mapped on the same contig
			if( align.IsFirstMate() && align.IsMateMapped() && align.RefID == align.MateRefID )
			{
				int32_t iSize;
				
				if( startRead < startMate )
				{
					iSize = (startMate + align.Length) - startRead;
					if( iSize < _minInsert[libId] || iSize > _maxInsert[libId] ) continue;
					
					// if the read and its mate are properly oriented update mean and std
					if( !align.IsReverseStrand() && align.IsMateReverseStrand() )
					{
						if(_isize_count[libId] == 1)
						{
							_isize_mean[libId] = iSize;
							_isize_std[libId] = 0;
							_isize_count[libId]++;
						}
						else
						{
							double oldMean = _isize_mean[libId];
							double oldStd = _isize_std[libId];
							
							_isize_mean[libId] = oldMean + (iSize - oldMean)/double(_isize_count[libId]);
							_isize_std[libId] = oldStd + (_isize_count[libId]-1)*(iSize - oldMean)*(iSize - oldMean)/double(_isize_count[libId]);
							_isize_count[libId]++;
						}
					}
				}
				else
				{
					iSize = (startRead + alignmentLength) - startMate;
					if( iSize < _minInsert[libId] || iSize > _maxInsert[libId] ) continue;
					
					// if the read and its mate are properly oriented update mean and std
					if( align.IsReverseStrand() && !align.IsMateReverseStrand() )
					{
						if(_isize_count[libId] == 1)
						{
							_isize_mean[libId] = iSize;
							_isize_std[libId] = 0;
							_isize_count[libId]++;
						}
						else
						{
							double oldMean = _isize_mean[libId];
							double oldStd = _isize_std[libId];
							
							_isize_mean[libId] = oldMean + (iSize - oldMean)/double(_isize_count[libId]);
							_isize_std[libId] = oldStd + (_isize_count[libId]-1)*(iSize - oldMean)*(iSize - oldMean)/double(_isize_count[libId]);
							_isize_count[libId]++;
						}
					}
				}
			}
		}
		
		// compute standard deviation
		this->_isize_std[libId] = sqrt( _isize_std[libId] / double(_isize_count[libId]) );
		
		// compute library's mean coverage
		this->_coverage[libId] = (this->_asm_size != 0) ? this->_reads_len[libId] / ((double)this->_asm_size) : 0.0;
	}
	
	this->Rewind();
	
	return true;
}


void MultiBamReader::writeStatsToFile( const std::string &filename ) const
{
	std::ofstream ofs( filename.c_str() );

	for( size_t i=0; i < _bam_readers.size(); i++ )
	{
		ofs << _bam_readers[i]->GetFilename() << std::endl;
		ofs << _isize_mean[i] << "\t" << _isize_std[i] << "\t" << _coverage[i] << std::endl; //"\t" << _isize_count[i] << std::endl;
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
			std::cerr << "The number of inserts statistics does not match the number of BAM files.\nStats file: " << filename << std::endl;
			return idx;
		}

		if( bamfile != _bam_readers[idx]->GetFilename() )
		{
			std::cerr << "Error loading bam statistics file (corresponding BAM file not found): " << bamfile << std::endl;
			return idx;
		}
		
		_isize_mean[idx] = _isize_std[idx] = _coverage[idx] = 0.0;

		getline( ifs, data );
		std::stringstream ss(data);
		ss >> _isize_mean[idx] >> _isize_std[idx] >> _coverage[idx]; // >> _isize_count[idx];

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