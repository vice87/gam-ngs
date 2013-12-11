/*
 *  This file is part of GAM-NGS.
 *  Copyright (c) 2011 by Riccardo Vicedomini <rvicedomini@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Simone Scalabrin <scalabrin@appliedgenomics.org>,
 *  Lars Arverstad <lars.arvestad@scilifelab.se>,
 *  Alberto Policriti <policriti@appliedgenomics.org>,
 *  Alberto Casagrande <casagrande@appliedgenomics.org>
 *
 *  GAM-NGS is an evolution of a previous work (GAM) done by Alberto Casagrande,
 *  Cristian Del Fabbro, Simone Scalabrin, and Alberto Policriti.
 *  In particular, GAM-NGS has been adapted to work on NGS data sets and it has
 *  been written using GAM's software as starting point. Thus, it shares part of
 *  GAM's source code.
 *
 *  GAM-NGS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GAM-NGS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GAM-NGS.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "UtilityFunctions.hpp"

char * getPathBaseName( char *path )
{
    char *ptr = strrchr(path, '/');
    return ptr ? ptr + 1 : (char *)path;
}

std::string getPathBaseName( const std::string path )
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


void loadBamFileNames( const std::string &input_file, std::vector< std::string > &names, std::vector<int32_t> &minInsert, std::vector<int32_t> &maxInsert )
{
	names.resize(0);
	minInsert.resize(0);
	maxInsert.resize(0);

	int32_t min_is, max_is;
	std::string line1, line2; // read line
	std::ifstream ifs( input_file.c_str() );

	while( ifs.good() )
	{
		min_is = 0;
		max_is = 0;

		getline(ifs,line1);

		if( line1 != "" )
		{
			names.push_back(line1);
			getline(ifs,line2);

			if( line2 != "" )
			{
				std::stringstream ss( line2 );
				ss >> min_is >> max_is;
				minInsert.push_back( min_is );
				maxInsert.push_back( max_is );
			}
			else
			{
				minInsert.push_back( min_is ); // min_is = 0
				maxInsert.push_back( max_is ); // max_is = 0
			}
		}
	}

	ifs.close();
}


void loadFileNames( const std::string &input_file, std::vector< std::string > &names )
{
	std::string line; // read line
	std::ifstream ifs( input_file.c_str() );

	while( ifs.good() )
	{
		getline(ifs,line);
		if( line != "" ) names.push_back(line);
	}

	ifs.close();
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


void mem_usage( double& vm_usage, double& resident_set )
{
	vm_usage     = 0.0;
	resident_set = 0.0;

	std::ifstream stat_stream( "/proc/self/stat", std::ios_base::in );

	// dummy vars for leading entries in stat that we don't care about
	std::string pid, comm, state, ppid, pgrp, session, tty_nr;
	std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	std::string utime, stime, cutime, cstime, priority, nice;
	std::string O, itrealvalue, starttime;

	// the two needed fields
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	>> utime >> stime >> cutime >> cstime >> priority >> nice
	>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}


void print_mem_usage()
{
	double vm_usage, rss_usage;
	mem_usage(vm_usage, rss_usage);

	std::string vm_suff = "KB", rss_suff = "KB";

	if( vm_usage > 1024 ) // possibly compute memory usage in MB or GB
	{
		vm_usage = vm_usage / 1024;
		if( vm_usage <= 1024 ) vm_suff = "MB";
		if( vm_usage > 1024 ){ vm_usage = vm_usage / 1024; vm_suff = "GB"; }
	}

	if( rss_usage > 1024 ) // possibly compute memory usage in MB or GB
	{
		rss_usage = rss_usage / 1024;
		if( rss_usage <= 1024 ) rss_suff = "MB";
		if( rss_usage > 1024 ){ rss_usage = rss_usage / 1024; rss_suff = "GB"; }
	}

	std::cerr << "[print_mem_usage] vm = " << vm_usage << vm_suff << "; rss = " << rss_usage << rss_suff << std::endl;
}


