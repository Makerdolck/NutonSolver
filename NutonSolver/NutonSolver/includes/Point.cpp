
#include "Point.h"

Point::Point(double _x, double _y, size_t _ID, bool _fixed)
{
	x		= _x;
	y		= _y;
	dx		= 0;
	dy		= 0;
	ID		= _ID;
	fixed	= _fixed;
	checked = false;
}

Point::~Point()
{
}