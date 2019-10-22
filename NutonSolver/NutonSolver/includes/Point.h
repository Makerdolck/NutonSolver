#pragma once

class Point
{
public:
	double	x,
			y;

	double	dx,
			dy;

	bool	fixed;

public:
	Point(double _x = 0, double _y = 0);
	~Point();

};
