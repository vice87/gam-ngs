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

char * getPathBaseName( char *path );

std::string getPathBaseName( std::string path );

std::string getBaseFileName( std::string filename );

std::string formatTime( time_t seconds );

int getMaxRSS(int64_t *maxrsskb);

#endif	/* UTILITYFUNCTIONS_HPP */

