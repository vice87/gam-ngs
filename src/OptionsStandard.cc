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

#include "OptionsStandard.hpp"
#include <sys/stat.h>

namespace options {

bool OptionsStandard::process(int argc, char *argv[])
{
	struct stat st;

	this->argc = argc;
	this->argv = argv;

	// PROCESS PARAMETERS
	std::stringstream ss;
	ss << "gam-ngs" << std::endl << std::endl << "Allowed options";

	po::options_description desc(ss.str().c_str());
	desc.add_options()
		// commands
		("help", "produce help message")
			//("version", "print version and exit")
		("create-blocks", "finds blocks given two alignments")
		("merge", "merge two assemblies, given a .blocks file")

		// input
		("master-bam", po::value< std::string >(), "coordinate sorted BAM file of the master assembly")
		("slave-bam", po::value< std::string >(), "coordinate sorted BAM file of the slave assembly")
		("master-isize", po::value< std::string >(), "insert size statistics file corresponding to master assembly")
		("slave-isize", po::value< std::string >(), "insert size statistics file corresponding to slave assembly")
                ("master-mp-bam", po::value< std::string >(), "coordinate sorted mate-pairs BAM file of the master assembly" )
                ("slave-mp-bam", po::value< std::string >(), "coordinate sorted mate-pairs BAM file of the slave assembly" )
				("master-mp-isize", po::value< std::string >(), "insert size statistics file corresponding to master assembly MP alignments")
				("slave-mp-isize", po::value< std::string >(), "insert size statistics file corresponding to slave assembly MP alignments")
		("master-namesorted-bam", po::value< std::string >(), "name sorted BAM file of the master assembly")
		("slave-namesorted-bam", po::value< std::string >(), "name sorted BAM file of the slave assembly")
		("blocks-file", po::value< std::string >(), "file .blocks created with the --create-blocks option")
		("master-fasta", po::value< std::string >(), "fasta file of the master assembly")
		("slave-fasta", po::value< std::string >(), "fasta file of the slave assembly")

                ("reads-prefix", po::value< std::string >(), "common prefix of all reads" )
		("min-block-size", po::value<int>(), "number of reads needed to build a block [default 50]")
		("threads", po::value<int>(), "number of threads [default 1]")
		("coverage-filter", po::value<double>(), "coverage filter threshold (default 0.75)")

		// output
		("output-prefix", po::value< std::string >(), "prefix of the output file")
		;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error error) {
		std::cerr <<  error.what() << std::endl;
		std::cerr << "Try \"--help\" for help" << std::endl;
		exit(2);
	}

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		exit(0);
	}

	/*if (vm.count("version")) {
		DEFAULT_CHANNEL << package_description() << endl;
		exit(0);
	}*/


	// COMMON PARAMETERS

	if( not(vm.count("create-blocks") or vm.count("merge")) )
	{
		std::cerr << "One between --create-blocks and --merge must be specified" << std::endl;
		std::cerr << "Try \"--help\" for help" << std::endl;
		exit(1);
	}

	if( vm.count("create-blocks") and vm.count("merge") )
	{
		std::cerr << "Only one between --create-blocks and --merge must be specified" << std::endl;
		exit(1);
	}

	if( not( vm.count("master-bam") and vm.count("slave-bam") ) )
	{
		std::cerr << "Either --master-bam and --slave-bam are required" << std::endl;
		exit(1);
	}

	if( not(vm.count("output-prefix")) )
	{
		std::cerr << "--output-prefix parameter is mandatory" << std::endl;
		exit(1);
	}

	masterBamFile = vm["master-bam"].as< std::string >();
	slaveBamFile = vm["slave-bam"].as< std::string >();

	// Check for master bam file existence */
	if( stat(masterBamFile.c_str(),&st) != 0 )
	{
		std::cerr << "Master BAM file " << masterBamFile << " doesn't exist." << std::endl;
		exit(1);
	}

	// Check for slave bam files existence */
	if( stat(slaveBamFile.c_str(),&st) != 0 )
	{
		std::cerr << "Slave BAM file " << slaveBamFile << " doesn't exist" << std::endl;
		exit(1);
	}

	if( vm.count("min-block-size") )
	{
		minBlockSize = vm["min-block-size"].as<int>();
		if( minBlockSize < 1 ) std::cerr << "warning: min-block-size is less than 1" << std::endl;
	}

	outputFilePrefix = vm["output-prefix"].as< std::string >();


	if( vm.count("create-blocks") ) // LOAD CREATE-BLOCKS PARAMETERS
	{
		program_mode = program_create_blocks;

		if( vm.count("master-namesorted-bam") )
		{
			masterNameSortedBamFile = vm["master-namesorted-bam"].as< std::string >();

			// Check for file existence */
			if( stat(masterNameSortedBamFile.c_str(),&st) != 0 )
			{
				std::cerr << "Master name-sorted BAM file " << masterNameSortedBamFile << " doesn't exist." << std::endl;
				exit(1);
			}
		}

		if( vm.count("slave-namesorted-bam") )
		{
			slaveNameSortedBamFile = vm["slave-namesorted-bam"].as< std::string >();

			// Check for file existence
			if( stat(slaveNameSortedBamFile.c_str(),&st) != 0 )
			{
				std::cerr << "Slave name-sorted BAM file " << slaveNameSortedBamFile << " doesn't exist." << std::endl;
				exit(1);
			}
		}

		if( vm.count("reads-prefix") )
		{
			this->readsPrefix = vm["reads-prefix"].as< std::string >();
		}
	}
	else if( vm.count("merge") ) // LOAD MERGE PARAMETERS
	{
		program_mode = program_merge;

		if( not( vm.count("master-fasta") and vm.count("slave-fasta") ) )
		{
			std::cerr << "Either --master-fasta and --slave-fasta are required" << std::endl;
			exit(1);
		}

		if( not( vm.count("blocks-file") ) )
		{
			std::cerr << "--blocks-file parameter is mandatory" << std::endl;
			exit(1);
		}

		masterFastaFile = vm["master-fasta"].as< std::string >();
		slaveFastaFile = vm["slave-fasta"].as< std::string >();

		// Check for master fasta file existence */
		if( stat(masterFastaFile.c_str(),&st) != 0 )
		{
			std::cerr << "Master FASTA file " << masterFastaFile << " doesn't exist." << std::endl;
			exit(1);
		}

		// Check for slave fasta files existence */
		if( stat( slaveFastaFile.c_str(),&st) != 0 )
		{
			std::cerr << "Slave BAM file " << slaveFastaFile << " doesn't exist" << std::endl;
			exit(1);
		}

		blocksFile = vm["blocks-file"].as< std::string >();

		// Check for blocks file existence */
		if( stat(blocksFile.c_str(),&st) != 0 )
		{
			std::cerr << "Blocks file " << blocksFile << " doesn't exist." << std::endl;
			exit(1);
		}

		if( vm.count("master-isize") ) masterISizeFile = vm["master-isize"].as< std::string >();
		if( vm.count("slave-isize") ) slaveISizeFile = vm["slave-isize"].as< std::string >();

                if( vm.count("master-mp-bam") ) masterMpBamFile = vm["master-mp-bam"].as< std::string >();
                if( vm.count("slave-mp-bam") ) slaveMpBamFile = vm["slave-mp-bam"].as< std::string >();
				if( vm.count("master-mp-isize") ) masterMpISizeFile = vm["master-mp-isize"].as< std::string >();
				if( vm.count("slave-mp-isize") ) slaveMpISizeFile = vm["slave-mp-isize"].as< std::string >();

                if( not vm.count("min-block-size") ) minBlockSize = 1;

		if( vm.count("threads") )
		{
			threadsNum = vm["threads"].as<int>();
			if( threadsNum < 1 ) threadsNum = 1;
		}

		if( vm.count("coverage-filter") )
		{
			if( coverageThreshold >= 0 ) coverageThreshold = vm["coverage-filter"].as<double>();
		}
	}

	return true;
}

} // end of namespace options
