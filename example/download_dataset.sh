#!/bin/bash

#  This file is part of GAM-NGS.
#  Copyright (c) 2011 by Riccardo Vicedomini <rvicedomini@appliedgenomics.org>,
#  Francesco Vezzi <vezzi@appliedgenomics.org>,
#  Simone Scalabrin <scalabrin@appliedgenomics.org>,
#  Lars Arverstad <lars.arvestad@scilifelab.se>,
#  Alberto Policriti <policriti@appliedgenomics.org>,
#  Alberto Casagrande <casagrande@appliedgenomics.org>
#
#  GAM-NGS is an evolution of a previous work (GAM) done by Alberto Casagrande,
#  Cristian Del Fabbro, Simone Scalabrin, and Alberto Policriti.
#  In particular, GAM-NGS has been adapted to work on NGS data sets and it has
#  been written using GAM's software as starting point. Thus, it shares part of
#  GAM's source code.
#
#  Moreover, GAM-NGS uses BamTools library to access BAM files.
#  BamTools's source code has been put in ./lib/bamtools-2.0.5/ folder,
#  in which its license can be found.
#
#  GAM-NGS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GAM-NGS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with GAM-NGS.  If not, see <http://www.gnu.org/licenses/>.

echo "Downloading assemblies."
(curl -# http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Assembly.tgz > ./Assembly.tgz) && tar -xzf ./Assembly.tgz
if [ $? -eq 0 ] ; then echo "Data set successfully retrieved." ; else (echo "There was a problem in retrieving the data set.";exit) ; fi

echo -e "\n"

echo "Downloading reads data set..."
(curl -# http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.allpathsCor.tgz > ./Data.allpathsCor.tgz) && tar -xzf ./Data.allpathsCor.tgz
if [ $? -eq 0 ] ; then echo "Data set successfully retrieved." ; else (echo "There was a problem in retrieving the data set.";exit) ; fi

