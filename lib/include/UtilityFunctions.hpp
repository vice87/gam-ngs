/*
 * File:   UtilityFunctions.hpp
 * Author: riccardo
 *
 * Created on 26 ottobre 2011, 17.58
 */

#ifndef UTILITYFUNCTIONS_HPP
#define	UTILITYFUNCTIONS_HPP

#include <stdint.h>
#include <string>
#include <vector>

char * getPathBaseName( char *path );

std::string getPathBaseName( const std::string path );

std::string getBaseFileName( std::string filename );

std::string formatTime( time_t seconds );

int getMaxRSS(int64_t *maxrsskb);

void loadFileNames( const std::string &input_file, std::vector< std::string > &names );

#endif	/* UTILITYFUNCTIONS_HPP */

