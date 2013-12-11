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

#include "OptionsCreate.hpp"
#include <sys/stat.h>

namespace options {

bool OptionsCreate::process(int argc, char *argv[])
{
	struct stat st;

	this->argc = argc;
	this->argv = argv;

	std::stringstream ss;
	ss << "\ngam-create: GAM-NGS executable for building blocks, given two alignments. Allowed options";

	po::options_description desc(ss.str().c_str());
	desc.add_options()
		// commands
		("help", "produce this help message\n")
		//("version", "print version and exit")

		// input
		("master-bam", po::value< std::string >(), "coordinate-sorted PE alignments of the master assembly")
		("slave-bam", po::value< std::string >(), "coordinate-sorted PE alignments of the slave assembly")

		//("master-namesorted-bam", po::value< std::string >(), "name sorted BAM file of the master assembly")
		//("slave-namesorted-bam", po::value< std::string >(), "name sorted BAM file of the slave assembly")

        ("min-block-size", po::value<int>(), "minimum number of reads needed to build a block (optional) [default=50]")
		//("threads", po::value<int>(), "number of threads [default 1]")

		// output
		("output", po::value< std::string >(), "output-file's prefix (optional) [default=out]")
		;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error error) {
		std::cerr <<  error.what() << std::endl;
		std::cerr << "Try \"--help\" for help." << std::endl;
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


	// INPUT PARAMETERS

	if( not( vm.count("master-bam") and vm.count("slave-bam") ) )
	{
		std::cerr << "Both --master-bam and --slave-bam options are required." << std::endl;
		exit(1);
	}

	masterBamFile = vm["master-bam"].as< std::string >();
	slaveBamFile = vm["slave-bam"].as< std::string >();

	// Check for master bam file existence */
	if( stat(masterBamFile.c_str(),&st) != 0 )
	{
		std::cerr << "Master BAM file " << masterBamFile << " does not exist." << std::endl;
		exit(1);
	}

	// Check for slave bam files existence */
	if( stat(slaveBamFile.c_str(),&st) != 0 )
	{
		std::cerr << "Slave BAM file " << slaveBamFile << " does not exist" << std::endl;
		exit(1);
	}

	if( vm.count("min-block-size") )
	{
		minBlockSize = vm["min-block-size"].as<int>();
		if( minBlockSize < 1 ) std::cerr << "WARNING: min-block-size is less than 1" << std::endl;
	}

	// OUTPUT
	if( vm.count("output") )
	{
		outputFilePrefix = vm["output"].as< std::string >();
	}

	return true;
}

} // end of namespace options
