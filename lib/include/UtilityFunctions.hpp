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

#ifndef UTILITYFUNCTIONS_HPP
#define UTILITYFUNCTIONS_HPP

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

char* getPathBaseName( char *path );

std::string getPathBaseName( const std::string path );

std::string getBaseFileName( std::string filename );

std::string formatTime( time_t seconds );

int getMaxRSS(int64_t *maxrsskb);

void mem_usage( double& vm_usage, double& resident_set );

void print_mem_usage();

void loadBamFileNames(
	const std::string &input_file,
	std::vector< std::string > &names,
	std::vector<int32_t> &minInsert,
	std::vector<int32_t> &maxInsert
);

void loadFileNames(
	const std::string &input_file,
	std::vector< std::string > &names
);

#endif	/* UTILITYFUNCTIONS_HPP */
