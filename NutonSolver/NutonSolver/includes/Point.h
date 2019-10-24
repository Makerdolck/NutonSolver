#pragma once

class Point
{
public:
	double	x,
			y;

	double	dx,
			dy;

	bool	fixed;
	bool 	checked;

public:
	Point(double _x = 0, double _y = 0, bool _fixed = false);
	~Point();

};
