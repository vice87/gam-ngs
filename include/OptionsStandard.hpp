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

#ifndef OPTIONS_STANDARD_H_
#define OPTIONS_STANDARD_H_

#include "Options.hpp"

namespace options {

class OptionsStandard : public Options{
public:
	OptionsStandard() : Options() { }
	OptionsStandard(int argc, char *argv[]) : Options() {
		if (not process(argc,argv))
			exit(2);
	}
	virtual ~OptionsStandard() { }

	bool process(int argc, char *argv[]);
};

} // end of namespace options

#endif /* OPTIONS_STANDARD_H_ */
