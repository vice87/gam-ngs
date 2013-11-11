#!/bin/bash

echo "Downloading assemblies."
(curl -# http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Assembly.tgz > ./Assembly.tgz) && tar -xzf ./Assembly.tgz
if [ $? -eq 0 ] ; then echo "Data set successfully retrieved." ; else (echo "There was a problem in retrieving the data set.";exit) ; fi

echo -e "\n"

echo "Downloading reads data set..."
(curl -# http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.allpathsCor.tgz > ./Data.allpathsCor.tgz) && tar -xzf ./Data.allpathsCor.tgz
if [ $? -eq 0 ] ; then echo "Data set successfully retrieved." ; else (echo "There was a problem in retrieving the data set.";exit) ; fi

