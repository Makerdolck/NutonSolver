
#include "Point.h"

Point::Point(double _x, double _y)
{
	x		= _x;
	y		= _y;
	dx		= 0;
	dy		= 0;
	fixed	= false;
}

Point::~Point()
{
}