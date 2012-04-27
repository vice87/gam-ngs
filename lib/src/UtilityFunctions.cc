
#include "UtilityFunctions.hpp"
#include <string.h>
#include <sstream>
#include <stdio.h>

char * getPathBaseName( char *path )
{
    char *ptr = strrchr(path, '/');
    return ptr ? ptr + 1 : (char *)path;
}

std::string getPathBaseName( std::string path )
{
    size_t found = path.rfind ( '/' );
    return (found != std::string::npos) ? path.substr(found+1) : std::string( path.c_str() );
}

std::string getBaseFileName( std::string filename )
{
    size_t found = filename.rfind('.');
    return filename.substr(0,found);
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

int getMaxRSS(int64_t *maxrsskb)
{
	int len = 0;
	int srtn = 0;
	char procf[257] = { "" };
	FILE *fp = NULL;
	char line[2001] = { "" };
	char crap[2001] = { "" };
	char units[2001] = { "" };
	int64_t maxrss = 0L;

	if(maxrsskb == NULL){
		return -1;
	}

	sprintf(procf,"/proc/%d/status",getpid());

	fp = fopen(procf, "r");
	if(fp == NULL){
		return -1;
	}

	while(fgets(line, 2000, fp) != NULL){
		if(strncasecmp(line,"VmPeak:",7) == 0){
			len = (int)strlen(line);
			line[len-1] = '\0';
	srtn = sscanf(line,"%s%ld%s",crap,&maxrss,units);
	if(srtn == 2){
		*maxrsskb = maxrss / 1024L;
	}else if(srtn == 3){
		if( (strcasecmp(units,"B") == 0) || (strcasecmp(units,"BYTES") == 0) ){
			*maxrsskb = maxrss / 1024L;
		}else if( (strcasecmp(units,"k") == 0) || (strcasecmp(units,"kB") == 0) ){
			*maxrsskb = maxrss * 1L;
		}else if( (strcasecmp(units,"m") == 0) || (strcasecmp(units,"mB") == 0) ){
			*maxrsskb = maxrss * 1024L;
		}else if( (strcasecmp(units,"g") == 0) || (strcasecmp(units,"gB") == 0) ){
			*maxrsskb = maxrss * 1024L * 1024L;
		}else{
			*maxrsskb = maxrss * 1L;
		}
	}
	break;
		}
	}

	fclose(fp);

	return 0;
}