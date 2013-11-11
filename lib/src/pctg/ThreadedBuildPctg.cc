/*
 * File:   ThreadedBuildPctg.code.hpp
 * Author: vice
 *
 * Created on 12 giugno 2011, 22.49
 */

#include <unistd.h>

#include "OptionsMerge.hpp"

#include "bam/MultiBamReader.hpp"
#include "graphs/AssemblyGraph.hpp"
#include "pctg/ThreadedBuildPctg.hpp"
#include "pctg/BuildPctgFunctions.hpp"

using namespace options;

extern OptionsMerge g_options;

extern MultiBamReader masterBam;
extern MultiBamReader masterMpBam;
extern MultiBamReader slaveBam;
extern MultiBamReader slaveMpBam;


CompactAssemblyGraph*
ThreadedBuildPctg::extractNextPctg()
{
	CompactAssemblyGraph *output = NULL;

    pthread_mutex_lock(&(this->_mutex));	

	while( _nextPctg < _graphs.size() )
	{
		if( boost::num_vertices(*(_graphs[_nextPctg])) == 0 )
		{
			_nextPctg++;
			continue;
		}
		else
		{
			output = _graphs[_nextPctg];
			_nextPctg++;
			break;
		}
	}

	pthread_mutex_unlock(&(this->_mutex));
	return output;
}

IdType
ThreadedBuildPctg::readPctgNumAndIncrease()
{
    IdType curPctgNum;

	pthread_mutex_lock(&(this->_mutexPctgNumInc));
    curPctgNum = this->_pctgNum;
    this->_pctgNum++;
	pthread_mutex_unlock(&(this->_mutexPctgNumInc));

    return curPctgNum;
}

void ThreadedBuildPctg::incProcBlocks( uint64_t num, uint64_t tid )
{
    pthread_mutex_lock(&(this->_mutexProcBlocks));

    this->_procBlocks += num;

    uint32_t perc = (100 * this->_procBlocks) / double(this->_totBlocks);
    if( this->_lastPerc < perc )
    {
		this->_lastPerc = perc;
		
		if( isatty(fileno(stdout)) )
		{
			std::cout << "\r[merge] Merging contigs " << perc << "% done." << std::flush;
			if( perc >= 100 ) std::cout << std::endl;
		}
		else
		{
			if( perc % 5 == 0 ) std::cout << "[merge] Merging contigs " << perc << "% done." << std::endl;
		}
    }

    pthread_mutex_unlock(&(this->_mutexProcBlocks));
}


ThreadedBuildPctg::ThreadedBuildPctg( 
	const std::list< CompactAssemblyGraph* > &graphsList, 
	const RefSequence &masterRef, 
	const RefSequence &slaveRef ) 
:
	_masterRef(masterRef), _slaveRef(slaveRef),
	_pctgNum(0), _nextPctg(0), _procBlocks(0), _totBlocks(0)
{
	(this->_graphs).resize( graphsList.size() );

    uint64_t b = 0;
    for( std::list< CompactAssemblyGraph* >::const_iterator it = graphsList.begin(); it != graphsList.end(); it++ )
    {
		this->_totBlocks += boost::num_vertices(**it);
		this->_graphs[b] = *it;
		
		b++;
    }

    pthread_mutex_init( &(this->_mutex), NULL );
    pthread_mutex_init( &(this->_mutexProcBlocks), NULL );
    pthread_mutex_init( &(this->_mutexPctgNumInc), NULL );

    pthread_mutex_init( &(this->_mutexMasterBam), NULL );
    pthread_mutex_init( &(this->_mutexSlaveBam), NULL );
}


std::list< PairedContig >*
ThreadedBuildPctg::run()
{
    this->_pctgsDone = 0;
    this->_lastPerc = 0;
	this->_procBlocks = 0;

	int threadsNum = (g_options.threadsNum > 0) ? g_options.threadsNum : 1;

	thread_arg_t* threads_argv[threadsNum];
	pthread_t threads[threadsNum];

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for( size_t i=0; i < threadsNum; i++ )
	{
		threads_argv[i] = new thread_arg_t;
		threads_argv[i]->tbp = this;
		threads_argv[i]->tid = i;

		threads_argv[i]->output = new std::list< PairedContig >;

		pthread_create( &threads[i], &attr, buildPctgThread, (void*)threads_argv[i] );
		std::cerr << "[build pctg] Thread " << i << " created." << std::endl;
	}

    pthread_attr_destroy(&attr);

    // wait on the other threads
    for( size_t i=0; i < threadsNum; i++ )
        pthread_join( threads[i], NULL );

	std::list< PairedContig > *outPctgList = new std::list< PairedContig >;

	// join PairedContig lists produced by each thread
	for( size_t i=0; i < threadsNum; i++ )
		outPctgList->splice( outPctgList->end(), *(threads_argv[i]->output) );

	// free dynamically allocated graphs
	for( size_t i=0; i < _graphs.size(); i++ ) 
	{
		if(_graphs[i] != NULL) delete _graphs[i];
	}

	// free dynamically allocated threads' arguments
	for( size_t i=0; i < threadsNum; i++ )
	{
		delete threads_argv[i]->output;
		delete threads_argv[i];
	}

	return outPctgList;
}


double ThreadedBuildPctg::computeZScore( MultiBamReader &multiBamReader, int32_t refID, uint32_t start, uint32_t end,	bool isMaster )
{
	uint32_t libs,idx;
	BamReader* bamReader;
	double lib_isize_mean, lib_isize_std;

	double z_score = 0;
	uint32_t minInsertNum = 5;

	libs = multiBamReader.size();

	uint64_t isize_num;

	/*for( size_t i=0; i < libs; i++ ) // find best library
	 {              *
	 if( i == 0 )
	 {
		 isize_num = multiBamReader.getISizeNum(i);
		 idx = i;
		 lib_isize_mean = multiBamReader.getISizeMean(i);
		 lib_isize_std = multiBamReader.getISizeStd(i);
		 bamReader = multiBamReader.getBamReader(i);
	 }
	 else if( multiBamReader.getISizeNum(i) > isize_num )
	 {
		 isize_num = multiBamReader.getISizeNum(i);
		 idx = i;
		 lib_isize_mean = multiBamReader.getISizeMean(i);
		 lib_isize_std = multiBamReader.getISizeStd(i);
		 bamReader = multiBamReader.getBamReader(i);
	 }
	 }*/

	if( libs == 0 ) return double(0);

	isize_num = multiBamReader.getISizeNum(0);
	idx = 0;
	lib_isize_mean = multiBamReader.getISizeMean(0);
	lib_isize_std = multiBamReader.getISizeStd(0);
	bamReader = multiBamReader.getBamReader(0);

	if( lib_isize_std == 0 ) return double(0);

	unsigned int times_std = 3;
	uint32_t min_insert = ( lib_isize_mean > times_std*lib_isize_std ) ? lib_isize_mean - times_std*lib_isize_std : 0;
	uint32_t max_insert = lib_isize_mean + times_std*lib_isize_std;

	bamReader->SetRegion( refID, start, refID, end+1 );

	BamAlignment align;
	uint64_t inserts=0, spanCov=0;
	int32_t nh, xt; // molteplicità delle read (nh->standard, xt->bwa)

	while( bamReader->GetNextAlignmentCore(align) ) // for each read in the region
	{
		if( !align.IsMapped() || align.Position < 0 || align.IsDuplicate() || !align.IsPrimaryAlignment() || align.IsFailedQC() ||
			!align.IsMateMapped() || align.RefID != align.MateRefID ) continue;

		int32_t read_start = align.Position;
		int32_t read_end = align.GetEndPosition() - 1;
		int32_t read_len = read_end - read_start + 1;
		int32_t mate_start = align.MatePosition;
		int32_t mate_end = align.MatePosition + read_len - 1;

		if( read_start < start || read_end > end ) continue;
		if( mate_start < start || mate_end > end ) continue;

		align.BuildCharData();
		if( !align.GetTag(std::string("NH"),nh) ) nh = 1;	// standard SAM format field
        if( !align.GetTag(std::string("XT"),xt) ) xt = 'U'; // bwa field
		bool is_uniq_mapped = (nh == 1 && xt == 'U');

		if( !is_uniq_mapped ) continue;

		if( align.IsFirstMate() )
		{
			if( read_start < mate_start )
			{
				int32_t i_size = (mate_start + read_len) - read_start;
				if( i_size < min_insert || i_size > max_insert ) continue;

				inserts++;
				spanCov += i_size;
			}
			else
			{
				int32_t i_size = read_end - mate_start + 1;
				if( i_size < min_insert || i_size > max_insert ) continue;

				inserts++;
				spanCov += i_size;
			}
		}
	}

	if( inserts > minInsertNum )
	{
		double localMean = spanCov/(double)inserts;
		z_score   = (localMean - lib_isize_mean)/(double)(lib_isize_std/sqrt(inserts));
	}

	return z_score;
}


void*
buildPctgThread(void* argv)
{
	// retrieve arguments of current thread
	
	typedef ThreadedBuildPctg::thread_arg_t thread_arg_t;

	thread_arg_t* thread_argv = (thread_arg_t*)argv;
    ThreadedBuildPctg *tbp = thread_argv->tbp;
	uint64_t tid = thread_argv->tid;

	std::list< PairedContig > *pctgList = thread_argv->output;
	CompactAssemblyGraph* cg = tbp->extractNextPctg();

	// process graphs
	while( cg != NULL )
	{
		try
        {
            buildPctg( tbp, *cg, tbp->_masterRef, tbp->_slaveRef, *pctgList );
        }
        catch(...) // this should not happen!
        {
            std::cerr << "Something unexpected happened processing graph " << cg->getId() << std::endl;
        }
		
		uint64_t cg_size = boost::num_vertices(*cg);
		tbp->incProcBlocks( cg_size, tid );
		
		cg->clear();
		cg = tbp->extractNextPctg();
	}

    pthread_exit((void *)0);
}
