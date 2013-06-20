/*
 * File:   UtilityFunctions.hpp
 * Author: riccardo
 *
 * Created on 26 ottobre 2011, 17.58
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
