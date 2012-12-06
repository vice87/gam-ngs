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

#include <iostream>
#include <iomanip>

#include "OptionsMerge.hpp"
#include "UtilityFunctions.hpp"
#include "Merge.hpp"

#include <boost/filesystem.hpp>

using namespace modules;

OptionsMerge g_options;

int main(int argc, char *argv[])
{
	if( not g_options.process(argc,argv) ) exit(2);

	Merge gamMerge;
	gamMerge.execute();
	
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

	std::cout << setiosflags(std::ios::fixed) << std::setprecision(2)
	          << "[gam-merge] MAX Memory used: " << maxrss << " " << maxrss_suff << std::endl;

    return 0;
}

