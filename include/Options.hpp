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

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <string>
#include <iostream>

namespace options {

class Options {
public:
	Options() { set_defaults(); }
	virtual ~Options() { }

	virtual bool process(int argc, char *argv[]) = 0;

	typedef enum
        {
            program_unknown,
            program_create_blocks,
            program_merge
        } __attribute__((packed)) program_mode_t;

	program_mode_t program_mode;

	int argc;
	char **argv;

	// input options
	std::string masterBamFile;
	std::string masterISizeFile;
	std::string slaveBamFile;
	std::string slaveISizeFile;

	std::string masterMpBamFile;
	std::string masterMpISizeFile;
	std::string slaveMpBamFile;
	std::string slaveMpISizeFile;

	std::string blocksFile;

	std::string masterFastaFile;
	std::string slaveFastaFile;

	int minBlockSize;
	int threadsNum;
	double coverageThreshold;
	bool noMultiplicityFilter;

	bool debug;

	bool outputGraphs;

	// unused input options
	std::string readsPrefix;
	std::string masterNameSortedBamFile;
	std::string slaveNameSortedBamFile;

	// output options
	std::string outputFilePrefix;

protected:
	void set_defaults();

};

} // end of namespace options

#endif /* OPTIONS_H_ */
