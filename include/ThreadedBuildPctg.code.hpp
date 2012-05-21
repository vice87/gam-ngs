/*
 * File:   ThreadedBuildPctg.code.hpp
 * Author: vice
 *
 * Created on 12 giugno 2011, 22.49
 */

#ifndef THREADEDBUILDPCTG_CODE_HPP
#define	THREADEDBUILDPCTG_CODE_HPP

#include "ThreadedBuildPctg.hpp"

#include "BuildPctgFunctions.code.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "graphs/GraphFunctions.hpp"

std::pair< std::vector<Block>, bool >
ThreadedBuildPctg::extractNextPctg()
{
    bool toDo;
    std::vector<Block> pctgVect;

    pthread_mutex_lock(&(this->_mutex));

    if( this->_blocksList.size() > 0 )
    {
        toDo = true;
        pctgVect = this->_blocksList.front();
        this->_blocksList.pop_front();

        this->_pctgsDone++;
        UIntType perc = 100 * ((RealType)(this->_pctgsDone)) / ((RealType)(this->_pctgsDone + this->_blocksList.size()));

        if( this->_lastPerc < perc )
        {
            this->_lastPerc = perc;
            if( perc % 5 == 0 ) std::cout << "[merging " << perc << "% done.]" << std::endl << std::flush;
        }
    }
    else
    {
        toDo = false;
    }

    pthread_mutex_unlock(&(this->_mutex));

    return std::pair< std::vector<Block>, bool >(pctgVect,toDo);
}

IdType
ThreadedBuildPctg::readPctgNumAndIncrease()
{
    IdType curPctgNum;

    pthread_mutex_lock(&(this->_mutex));
    curPctgNum = this->_pctgNum;
    this->_pctgNum++;
    pthread_mutex_unlock(&(this->_mutex));

    return curPctgNum;
}


void
ThreadedBuildPctg::memorizeNewPctg( std::list< PairedContig > &pctgList )
{
    pthread_mutex_lock(&(this->_mutexRemoveCtgId));

    //(this->_pctgPool)->set( pctg.name(), pctg ); // aggiungi il paired contig al rispettivo pool

    // segna i contigs del master assembly che sono stati inseriti in un paired contig
    std::list< PairedContig >::iterator pctg;
    for( pctg = pctgList.begin(); pctg != pctgList.end(); pctg++ )
    {
        std::map< std::pair<IdType,IdType>, ContigInPctgInfo >::const_iterator ctg_id;
        for( ctg_id = (pctg->getMasterCtgMap()).begin(); ctg_id != (pctg->getMasterCtgMap()).end(); ctg_id++ )
        {
            this->_removedCtgs[ (ctg_id->first).second ] = true;
        }
    }

    // aggiorna la lista dei paired contig trovati, inserendo quelli nuovi
    (this->_pctgList).splice( (this->_pctgList).end(), pctgList );

    pthread_mutex_unlock(&(this->_mutexRemoveCtgId));
}


ThreadedBuildPctg::ThreadedBuildPctg(
        const std::list<std::vector<Block> >& blocksList,
        HashContigMemPool* pctgPool,
        const ExtContigMemPool *masterPool,
        const ExtContigMemPool *slavePool,
		MultiBamReader &masterBam,
		std::vector< MultiBamReader > &slaveBams,
        const BamTools::RefVector *masterRefVector,
        const std::vector< BamTools::RefVector > *slaveRefVector,
		const Options &options ):
                _pctgNum(0),
                _blocksList(blocksList),
                _pctgPool(pctgPool),
                _masterPool(masterPool),
                _slavePool(slavePool),
                _masterBam(masterBam),
                _slaveBams(slaveBams),
                _masterRefVector(masterRefVector),
                _slaveRefVector(slaveRefVector),
                _removedCtgs(masterRefVector->size(),false),
                _options(options)
{
    //this->_pctgNum = 0;
    //this->_pctgList = std::list< PairedContig >();
    //this->_removedCtgs = std::vector<bool>( masterRefVector->size(), false );

    this->_pctgPool->clear();

	pthread_mutex_init( &(this->_mutexRemoveCtgId), NULL );
    pthread_mutex_init( &(this->_mutex), NULL );

	pthread_mutex_init( &(this->_mutexMasterBam), NULL );

	_mutexSlaveBams.resize( _slaveBams.size() );
	for( size_t i=0; i < _mutexSlaveBams.size(); i++ ) pthread_mutex_init( &(this->_mutexSlaveBams[i]), NULL );
}


std::pair< std::list<PairedContig>, std::vector<bool> >
ThreadedBuildPctg::run()
{
    this->_pctgsDone = 0;
    this->_lastPerc = 0;

	int threadsNum = _options.threadsNum;

	pthread_t threads[threadsNum];

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for( size_t i=0; i < threadsNum; i++ )
        pthread_create( &threads[i], &attr, buildPctgThread, (void*)this );

    pthread_attr_destroy(&attr);

    // wait on the other threads
    for( size_t i=0; i < threadsNum; i++ )
        pthread_join( threads[i], NULL );

    return std::pair< std::list<PairedContig>, std::vector<bool> >( this->_pctgList, this->_removedCtgs );
}


double ThreadedBuildPctg::computeZScore( const std::pair<uint64_t,uint64_t> &ctgId, uint32_t start, uint32_t end,	bool isMaster )
{
	uint32_t libs,idx;
	BamReader* bamReader;
	double lib_isize_mean, lib_isize_std;

	double Z_stats = 0;
	uint32_t minInsertNum = 5;

	if( isMaster )
	{
		libs = _masterBam.size();

		uint64_t isize_num;
		for( size_t i=0; i < libs; i++ ) // find best library
		{
			if( i == 0 )
			{
				isize_num = _masterBam.getISizeNum(i);
				idx = i;
				lib_isize_mean = _masterBam.getISizeMean(i);
				lib_isize_std = _masterBam.getISizeStd(i);
				bamReader = _masterBam.getBamReader(i);
			}
			else if( _masterBam.getISizeNum(i) > isize_num )
			{
				isize_num = _masterBam.getISizeNum(i);
				idx = i;
				lib_isize_mean = _masterBam.getISizeMean(i);
				lib_isize_std = _masterBam.getISizeStd(i);
				bamReader = _masterBam.getBamReader(i);
			}
		}
	}
	else // not isMaster
	{
		libs = _slaveBams.at(ctgId.first).size();

		uint64_t isize_num;
		for( size_t i=0; i < libs; i++ ) // find best library
		{
			if( i == 0 )
			{
				isize_num = _slaveBams.at(ctgId.first).getISizeNum(i);
				idx = i;
				lib_isize_mean = _slaveBams.at(ctgId.first).getISizeMean(i);
				lib_isize_std = _slaveBams.at(ctgId.first).getISizeStd(i);
				bamReader = _slaveBams.at(ctgId.first).getBamReader(i);
			}
			else if( _slaveBams.at(ctgId.first).getISizeNum(i) > isize_num )
			{
				isize_num = _slaveBams.at(ctgId.first).getISizeNum(i);
				idx = i;
				lib_isize_mean = _slaveBams.at(ctgId.first).getISizeMean(i);
				lib_isize_std = _slaveBams.at(ctgId.first).getISizeStd(i);
				bamReader = _slaveBams.at(ctgId.first).getBamReader(i);
			}
		}
	}

	// lock mutex for accessing the bam file
	if( isMaster ) pthread_mutex_lock(&(this->_mutexMasterBam)); else pthread_mutex_lock(&(this->_mutexSlaveBams[ctgId.first]));

	bamReader->SetRegion( ctgId.second, start, ctgId.second, end+1 );

	BamAlignment align;
	uint64_t inserts=0, spanCov=0;

	while( bamReader->GetNextAlignmentCore(align) ) // for each read in the region
	{
		if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ||
			!align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

		if( align.Position < start ) continue;

		if( (align.MatePosition < start && align.MatePosition > end) || align.IsFirstMate()  )
		{
			if( align.Position < align.MatePosition /*&& !align.IsReverseStrand() && align.IsMateReverseStrand()*/ )
			{
				uint32_t iSize = (align.MatePosition + align.Length) - align.Position;
				if( iSize < MIN_ISIZE || iSize > MAX_ISIZE ) continue;

				inserts++;
				spanCov += iSize;
			}
			else if( align.Position >= align.MatePosition /*&& align.IsReverseStrand() && !align.IsMateReverseStrand()*/ )
			{
				int32_t alignmentLength = align.GetEndPosition() - align.Position;
				uint32_t iSize = (align.Position + alignmentLength) - align.MatePosition;
				if( iSize < MIN_ISIZE || iSize > MAX_ISIZE ) continue;

				inserts++;
				spanCov += iSize;
			}
		}
	}

	// unlock mutex used for accessing the bam file
	if( isMaster ) pthread_mutex_unlock(&(this->_mutexMasterBam)); else pthread_mutex_unlock(&(this->_mutexSlaveBams[ctgId.first]));

	if( inserts > minInsertNum )
	{
		double localMean = spanCov/(double)inserts;
		Z_stats   = (localMean - lib_isize_mean)/(double)(lib_isize_std/sqrt(inserts));
	}

	return Z_stats;
}


double ThreadedBuildPctg::computeReadCoverage( std::pair<uint64_t,uint64_t> &ctgId, uint32_t start, uint32_t end, bool isMaster )
{
	BamReader* bamReader;
	double coverage = 0;

	if( isMaster ) bamReader = _masterBam.getBamReader(0); else bamReader = _slaveBams.at(ctgId.first).getBamReader(0);

	// lock mutex for accessing the bam file
	if( isMaster ) pthread_mutex_lock(&(this->_mutexMasterBam)); else pthread_mutex_lock(&(this->_mutexSlaveBams[ctgId.first]));

	bamReader->SetRegion( ctgId.second, start, ctgId.second, end+1 );

	BamAlignment align;
	int64_t lengthSum = 0;

	while( bamReader->GetNextAlignmentCore(align) ) // for each read in the region
	{
		if( !align.IsMapped() || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ) continue;

		// TODO: tenere conto solo della porzione di read che si sovrappone alla regione definita
		int32_t alignmentLength = align.GetEndPosition() - align.Position;
		lengthSum += alignmentLength;
	}

	// unlock mutex used for accessing the bam file
	if( isMaster ) pthread_mutex_unlock(&(this->_mutexMasterBam)); else pthread_mutex_unlock(&(this->_mutexSlaveBams[ctgId.first]));

	uint64_t region_len = end >= start ? end-start+1 : 0;

	if( region_len != 0 ) coverage = double(lengthSum)/double(region_len);

	return coverage;
}


void*
buildPctgThread(void* argv)
{
    //IdType pctgId = 0;
    //bool newPctgIdRequired = true;

    ThreadedBuildPctg *tbp( (ThreadedBuildPctg *)argv );
    std::pair< std::vector<Block>, bool > nextPctg = tbp->extractNextPctg();

    while( nextPctg.second )
    {
        AssemblyGraph graph( nextPctg.first );

        if( is_linear(graph) )
        {
            try
            {
                std::list< PairedContig > pctgList;

                buildPctg( tbp, graph, tbp->_masterPool, tbp->_slavePool,
                        tbp->_masterRefVector, tbp->_slaveRefVector, pctgList, tbp->_options );

                std::list< PairedContig >::iterator pctg = pctgList.begin();
                while( pctg != pctgList.end() )
                {
                    if( pctg->size() != 0 )
                    {
                        pctg->setId( tbp->readPctgNumAndIncrease() );
                        ++pctg;
                    }
                    else
                    {
                        pctg = pctgList.erase(pctg);
                    }
                }

                tbp->memorizeNewPctg( pctgList );

                //newPctgIdRequired = true;
                //std::cout << std::endl << std::flush;
            }
            catch(...)
            {
                std::cerr << "I cannot build a pctg from a cycle graph." << std::endl;
            }
        }

        nextPctg = tbp->extractNextPctg();
        graph.clear();
    }

    pthread_exit((void *)0);
}

#endif	/* THREADEDBUILDPCTG_CODE_HPP */

