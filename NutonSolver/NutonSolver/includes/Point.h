#pragma once

#include <cstring>


class Point
{
public:
	double	x,
			y;

	double	dx,
			dy;

	size_t	ID;

	bool	fixed;
	bool 	checked;

public:
	Point(double _x = 0, double _y = 0, size_t _ID = 0, bool _fixed = false);
	~Point();

};
