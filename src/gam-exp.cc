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

#include "OptionsStandard.hpp"
#include "UtilityFunctions.hpp"

#include "Merge.hpp"
#include "CreateBlocks.hpp"

using namespace modules;

int main(int argc, char *argv[])
{
    OptionsStandard options(argc,argv);

    if( options.program_mode == Options::program_create_blocks )
    {
        CreateBlocks createBlocks;
        createBlocks.execute(options);
    }
    else if( options.program_mode == Options::program_merge )
    {
        Merge gamMerge;
        gamMerge.execute(options);
    }
    else
    {
        std::cerr << "Program mode not found!" << std::endl;
        exit(2);
    }

    int64_t maxrsskb = 0L;
    if( getMaxRSS( &maxrsskb ) == 0 ) std::cout << "MAX Memory used: " << maxrsskb << std::endl;

    return 0;
}

