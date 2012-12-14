/*
 *  This file is part of rNA.
 *  Copyright (c) 2011 by Cristian Del Fabbro <delfabbro@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Alexandru Tomescu <alexandru.tomescu@uniud.it>, and
 *  Alberto Policriti <policriti@uniud.it>
 *
 *   rNA is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rNA is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with rNA.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "OptionsMerge.hpp"
#include <sys/stat.h>

namespace options {

bool OptionsMerge::process(int argc, char *argv[])
{
	struct stat st;

	this->argc = argc;
	this->argv = argv;

	// PROCESS PARAMETERS
	std::stringstream ss;
	ss << "\ngam-merge: GAM-NGS executable for merging two assemblies. Allowed options";

	po::options_description desc(ss.str().c_str());
	desc.add_options()
		// commands
		("help", "produce this help message\n")
		//("version", "print version and exit")
		
		// input
		("master-bam", po::value< std::string >(), "coordinate-sorted PE alignments of the master assembly")
		("slave-bam", po::value< std::string >(), "coordinate-sorted PE alignments of the slave assembly")
		//("master-isize", po::value< std::string >(), "insert size statistics file corresponding to master assembly")
		//("slave-isize", po::value< std::string >(), "insert size statistics file corresponding to slave assembly")
		("master-mp-bam", po::value< std::string >(), "coordinate sorted MP alignments of the master assembly (optional)" )
		("slave-mp-bam", po::value< std::string >(), "coordinate sorted MP alignments of the slave assembly (optional)" )
		//("master-mp-isize", po::value< std::string >(), "insert size statistics file corresponding to master assembly MP alignments")
		//("slave-mp-isize", po::value< std::string >(), "insert size statistics file corresponding to slave assembly MP alignments")
		//("master-namesorted-bam", po::value< std::string >(), "name sorted BAM file of the master assembly")
		//("slave-namesorted-bam", po::value< std::string >(), "name sorted BAM file of the slave assembly")
		("blocks-file", po::value< std::string >(), ".blocks file created with gam-create command")
		("master-fasta", po::value< std::string >(), "fasta file of the master assembly")
		("slave-fasta", po::value< std::string >(), "fasta file of the slave assembly")
		//("reads-prefix", po::value< std::string >(), "common prefix of all reads" )
		("min-block-size", po::value<int>(), "minimum number of reads of blocks to be loaded (optional) [default=5]")
		("threads", po::value<int>(), "number of threads (optional) [default=1]")
		("coverage-filter", po::value<double>(), "coverage filter threshold (optional) [default=0.75]")
		
		("output-graphs", "output graphs in gam_graphs sub-folder (debug)")

		// output
		("output", po::value< std::string >(), "output-files' prefix (optional) [default=out]")
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

	if (vm.count("help")) 
	{
		std::cout << desc << std::endl;
		std::cout << "Updated sources and documentation can be found at http://github.com/vice87/gam-ngs\n" << std::endl;
		exit(0);
	}

	/*if (vm.count("version")) {
		DEFAULT_CHANNEL << package_description() << endl;
		exit(0);
	}*/


	// both master/slave alignments have to be provided
	if( not( vm.count("master-bam") and vm.count("slave-bam") ) )
	{
		std::cerr << "Both --master-bam and --slave-bam parameters are mandatory." << std::endl;
		std::cerr << "Try \"--help\" for help" << std::endl;
		exit(1);
	}
	else
	{
		masterBamFile = vm["master-bam"].as< std::string >();
		slaveBamFile = vm["slave-bam"].as< std::string >();
		
		masterISizeFile = masterBamFile + ".isize";
		slaveISizeFile = slaveBamFile + ".isize";
		
		// Check for master bam file existence */
		if( stat(masterBamFile.c_str(),&st) != 0 )
		{
			std::cerr << "Master's PE-alignments file " << masterBamFile << " does not exist." << std::endl;
			exit(1);
		}
		
		// Check for slave bam files existence */
		if( stat(slaveBamFile.c_str(),&st) != 0 )
		{
			std::cerr << "Slave's PE-aligments file " << slaveBamFile << " does not exist." << std::endl;
			exit(1);
		}
	}
	
	if( vm.count("master-mp-bam") || vm.count("slave-mp-bam") ) 
	{
		// both the parameters have to be given, or none of them
		if( not( vm.count("master-mp-bam") and vm.count("slave-mp-bam") ) )
		{
			std::cerr << "Both --master-mp-bam and --slave-mp-bam have to be specified, or none of them." << std::endl;
			std::cerr << "Try \"--help\" for help" << std::endl;
			exit(1);
		}
		
		masterMpBamFile = vm["master-mp-bam"].as< std::string >();
		slaveMpBamFile = vm["slave-mp-bam"].as< std::string >();
		
		// Check for master MP alignments file existence */
		if( stat(masterMpBamFile.c_str(),&st) != 0 )
		{
			std::cerr << "Master's MP-alignments file " << masterMpBamFile << " does not exist." << std::endl;
			exit(1);
		}
		
		// Check for slave MP alignments file existence */
		if( stat(slaveMpBamFile.c_str(),&st) != 0 )
		{
			std::cerr << "Slave's MP-aligments file " << slaveMpBamFile << " does not exist." << std::endl;
			exit(1);
		}
		
		if( masterMpBamFile != "" ) masterMpISizeFile = masterMpBamFile + ".isize";
		if( slaveMpBamFile != "" ) slaveMpISizeFile = slaveMpBamFile + ".isize";
	}
	
	
	if( not( vm.count("blocks-file") ) )
	{
		std::cerr << "--blocks-file parameter is mandatory." << std::endl;
		std::cerr << "Try \"--help\" for help" << std::endl;
		exit(1);
	}
	else
	{
		blocksFile = vm["blocks-file"].as< std::string >();
		
		// Check for blocks file existence */
		if( stat(blocksFile.c_str(),&st) != 0 )
		{
			std::cerr << "Blocks' file " << blocksFile << " does not exist." << std::endl;
			exit(1);
		}
	}
	
	
	if( not( vm.count("master-fasta") and vm.count("slave-fasta") ) )
	{
		std::cerr << "Both --master-fasta and --slave-fasta parameters are mandatory." << std::endl;
		std::cerr << "Try \"--help\" for help" << std::endl;
		exit(1);
	}
	else
	{
		masterFastaFile = vm["master-fasta"].as< std::string >();
		slaveFastaFile = vm["slave-fasta"].as< std::string >();
		
		// Check for master fasta file existence //
		if( stat(masterFastaFile.c_str(),&st) != 0 )
		{
			std::cerr << "Master-assembly's fasta file " << masterFastaFile << " does not exist." << std::endl;
			exit(1);
		}
		
		// Check for slave fasta files existence //
		if( stat( slaveFastaFile.c_str(),&st) != 0 )
		{
			std::cerr << "Slave-assembly's fasta file " << slaveFastaFile << " does not exist" << std::endl;
			exit(1);
		}
	}
	
	
	if( vm.count("min-block-size") )
	{
		minBlockSize = vm["min-block-size"].as<int>();
		if( minBlockSize < 1 ) std::cerr << "warning: min-block-size is less than 1" << std::endl;
	}
	else
	{
		minBlockSize = 5;
	}
	
	
	if( vm.count("threads") )
	{
		threadsNum = vm["threads"].as<int>();
		if( threadsNum < 1 ) threadsNum = 1;
	}
	
	
	if( vm.count("coverage-filter") )
	{
		if( coverageThreshold >= 0 ) coverageThreshold = vm["coverage-filter"].as<double>();
	}
	
	
	if( vm.count("output-graphs") )
	{
		outputGraphs = true;
	}
	
	
	if( vm.count("output") )
	{
		outputFilePrefix = vm["output"].as< std::string >();
	}

	
	return true;
}

} // end of namespace options
