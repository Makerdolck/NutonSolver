
#include "Point.h"

Point::Point(double _x, double _y, bool _fixed)
{
	x		= _x;
	y		= _y;
	dx		= 0;
	dy		= 0;
	fixed	= _fixed;
	checked = false;
}

Point::~Point()
{
}