/* 
 * File:   UtilityFunctions.hpp
 * Author: riccardo
 *
 * Created on 26 ottobre 2011, 17.58
 */

#ifndef UTILITYFUNCTIONS_HPP
#define	UTILITYFUNCTIONS_HPP

#include <string>

char * getPathBaseName( char *path )
{
    char *ptr = strrchr(path, '/');
    return ptr ? ptr + 1 : (char *)path;
}

std::string getPathBaseName( std::string path )
{
    size_t found = path.rfind ( '/' );
    return (found != std::string::npos) ? path.substr(found) : std::string( path.c_str() );
}

std::string formatTime( time_t seconds )
{
    std::stringstream out;
    int h = seconds / 3600;
    int m = (seconds % 3600) / 60;
    int s = (seconds % 3600) % 60;
        
    if( h > 0 ) out << h << "h"; 
    if( m > 0 ) out << m << "m";
    out << s << "s";
    
    return out.str();
}

#endif	/* UTILITYFUNCTIONS_HPP */

