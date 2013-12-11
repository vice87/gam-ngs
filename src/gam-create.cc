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

#include <iostream>
#include <iomanip>

#include "OptionsCreate.hpp"
#include "UtilityFunctions.hpp"
#include "CreateBlocks.hpp"

#include <boost/filesystem.hpp>

using namespace modules;

OptionsCreate g_options;

int main(int argc, char *argv[])
{
	if( not g_options.process(argc,argv) ) exit(2);

	CreateBlocks createBlocks;
    createBlocks.execute();

    int64_t maxrsskb = 0L;
	getMaxRSS( &maxrsskb );

	double maxrss = maxrsskb;
	std::string maxrss_suff = "KB";

	if( maxrss > 1024 )
	{
		maxrss = maxrss / 1024;
		if( maxrss <= 1024 ) maxrss_suff = "MB";
		if( maxrss > 1024 ){ maxrss = maxrss / 1024; maxrss_suff = "GB"; }
	}

	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2)
	          << "[gam-create] MAX Memory used: " << maxrss << " " << maxrss_suff << std::endl;

	//print_mem_usage();

    return 0;
}

