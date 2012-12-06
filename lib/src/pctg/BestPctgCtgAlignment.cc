#include "pctg/BestPctgCtgAlignment.hpp"

BestPctgCtgAlignment::BestPctgCtgAlignment() {}


BestPctgCtgAlignment::BestPctgCtgAlignment(const MyAlignment& main, const bool isCtgReverse):
        _main(1,main),
        _left(100),
        _right(100),
        _isCtgReverse(isCtgReverse),
        _left_rev(false),
        _right_rev(false)
{}


BestPctgCtgAlignment::BestPctgCtgAlignment(const std::vector<MyAlignment>& main, const bool isCtgReverse):
	_main(main),
	_left(100),
	_right(100),
	_isCtgReverse(isCtgReverse),
	_left_rev(false),
	_right_rev(false)
{}


BestPctgCtgAlignment::BestPctgCtgAlignment(const BestPctgCtgAlignment& orig):
        _main(orig._main),
        _left(orig._left),
        _right(orig._right),
        _isCtgReverse(orig._isCtgReverse),
        _left_rev(orig._left_rev),
        _right_rev(orig._right_rev)
{}


BestPctgCtgAlignment::BestPctgCtgAlignment(
	const MyAlignment& main,
	const bool isCtgReverse,
	const MyAlignment& left,
	const MyAlignment& right,
	const bool left_rev,
	const bool right_rev
) :
	_main(1,main),
	_left(left),
	_right(right),
	_isCtgReverse(isCtgReverse),
	_left_rev(left_rev),
	_right_rev(right_rev)
{}


BestPctgCtgAlignment::BestPctgCtgAlignment(
	const std::vector<MyAlignment>& main,
	const bool isCtgReverse,
	const MyAlignment& left,
	const MyAlignment& right,
	const bool left_rev,
	const bool right_rev
) :
	_main(main),
	_left(left),
	_right(right),
	_isCtgReverse(isCtgReverse),
	_left_rev(left_rev),
	_right_rev(right_rev)
{}


const std::vector<MyAlignment>&
BestPctgCtgAlignment::main() const
{
    return this->_main;
}

const MyAlignment&
BestPctgCtgAlignment::left() const
{
	return this->_left;
}

const MyAlignment&
BestPctgCtgAlignment::right() const
{
	return this->_right;
}


double BestPctgCtgAlignment::main_homology() const
{
	double min_homology = 0;

	for( size_t i=0; i < this->_main.size(); i++ )
	{
		if( i==0 )
		{
			min_homology = _main[i].homology();
		}
		else if( _main[i].homology() < min_homology )
		{
			min_homology = _main[i].homology();
		}
	}

	return min_homology;
}


bool
BestPctgCtgAlignment::is_left_rev() const
{
	return this->_left_rev;
}

bool
BestPctgCtgAlignment::is_right_rev() const
{
	return this->_right_rev;
}


bool
BestPctgCtgAlignment::isCtgReversed() const
{
    return this->_isCtgReverse;
}

const BestPctgCtgAlignment&
BestPctgCtgAlignment::operator =(const BestPctgCtgAlignment& orig)
{
    this->_main = (orig._main);
	this->_left = (orig._left);
	this->_right = (orig._right);

    this->_isCtgReverse = (orig._isCtgReverse);
	this->_left_rev = (orig._left_rev);
	this->_right_rev = (orig._right_rev);

    return *this;
}


const MyAlignment&
BestPctgCtgAlignment::operator[]( const size_t &index ) const
{
	return this->_main[index];
}


uint64_t
BestPctgCtgAlignment::size() const
{
	return this->_main.size();
}