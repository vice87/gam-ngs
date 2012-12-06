/*
 * File:   n50.cc
 * Author: riccardo
 *
 * Created on 1 agosto 2011, 22.39
 */

#include <cstdlib>
#include <stdint.h>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <map>

std::string getBaseName( const char *p )
{
	std::string path = p;

	size_t found = path.rfind('/');
	if( found != std::string::npos ) path = path.substr( found+1 );

	return path;
}

int main(int argc, char** argv)
{
	if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <id1.txt> [<id2.txt> ...]" << std::endl;
        return 1;
    }

    uint64_t duplicates = 0;

	typedef std::multimap< std::string, uint32_t > seq_map_t;

    seq_map_t seq_map;
	std::ifstream ifs;

	std::cout << "[main] loading sequences identifiers" << std::endl;

    for( size_t i=1; i < argc; i++ )
	{
		ifs.open( argv[i] , std::ifstream::in );

		while( ifs.good() )
		{
			std::string id;
			getline(ifs,id);

			if( id.length() == 0 ) continue;

			seq_map.insert( std::pair<std::string,uint32_t>(id,i) );
		}

		ifs.close();
	}

	std::set< std::string > keys;
	std::pair< std::string, uint32_t > last;
	bool last_dup = false;

	std::cout << "[main] building keys set" << std::endl;

	seq_map_t::const_iterator it = seq_map.begin();
	while( it != seq_map.end() ){ keys.insert( it->first ); ++it; }

	std::cout << "[main] searching for duplicates\n" << std::endl;

	std::set< std::string >::const_iterator it2 = keys.begin();
	while( it2 != keys.end() )
	{
		seq_map_t::iterator it;
		std::pair< seq_map_t::iterator, seq_map_t::iterator > ret;

		if( seq_map.count(*it2) > 1 )
		{
			ret = seq_map.equal_range( *it2 );
			for( it=ret.first; it!=ret.second; ++it )
			{
				if( it == ret.first ) std::cerr << it->first << " duplicate in following files:" << std::endl;
				std::cerr << "\t" << getBaseName( argv[it->second] ) << std::endl;
			}
		}

		++it2;
	}

    return 0;
}

