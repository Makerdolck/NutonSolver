#pragma once

#include "Point.h"

class Constraint
{
private:
	Point	tmp_pointA1,
			tmp_pointA2,
			tmp_pointB1,
			tmp_pointB2;

public:
	Point	*pointA1,
			*pointA2,
			*pointB1,
			*pointB2;

	double	value;

	double (*func)(Point, Point, Point, Point, double);

	double (*func_dA1x)(Point, Point, Point, Point, double);
	double (*func_dA1y)(Point, Point, Point, Point, double);
	double (*func_dA2x)(Point, Point, Point, Point, double);
	double (*func_dA2y)(Point, Point, Point, Point, double);
	double (*func_dB1x)(Point, Point, Point, Point, double);
	double (*func_dB1y)(Point, Point, Point, Point, double);
	double (*func_dB2x)(Point, Point, Point, Point, double);
	double (*func_dB2y)(Point, Point, Point, Point, double);

	//	A1x_
	double (*func_dA1x_dA1x)(Point, Point, Point, Point, double);
	double (*func_dA1x_dA1y)(Point, Point, Point, Point, double);
	double (*func_dA1x_dA2x)(Point, Point, Point, Point, double);
	double (*func_dA1x_dA2y)(Point, Point, Point, Point, double);
	double (*func_dA1x_dB1x)(Point, Point, Point, Point, double);
	double (*func_dA1x_dB1y)(Point, Point, Point, Point, double);
	double (*func_dA1x_dB2x)(Point, Point, Point, Point, double);
	double (*func_dA1x_dB2y)(Point, Point, Point, Point, double);
	// A1y_
	double (*func_dA1y_dA1x)(Point, Point, Point, Point, double);
	double (*func_dA1y_dA1y)(Point, Point, Point, Point, double);
	double (*func_dA1y_dA2x)(Point, Point, Point, Point, double);
	double (*func_dA1y_dA2y)(Point, Point, Point, Point, double);
	double (*func_dA1y_dB1x)(Point, Point, Point, Point, double);
	double (*func_dA1y_dB1y)(Point, Point, Point, Point, double);
	double (*func_dA1y_dB2x)(Point, Point, Point, Point, double);
	double (*func_dA1y_dB2y)(Point, Point, Point, Point, double);
	//
	//	A2x_
	double (*func_dA2x_dA1x)(Point, Point, Point, Point, double);
	double (*func_dA2x_dA1y)(Point, Point, Point, Point, double);
	double (*func_dA2x_dA2x)(Point, Point, Point, Point, double);
	double (*func_dA2x_dA2y)(Point, Point, Point, Point, double);
	double (*func_dA2x_dB1x)(Point, Point, Point, Point, double);
	double (*func_dA2x_dB1y)(Point, Point, Point, Point, double);
	double (*func_dA2x_dB2x)(Point, Point, Point, Point, double);
	double (*func_dA2x_dB2y)(Point, Point, Point, Point, double);
	// A2y_
	double (*func_dA2y_dA1x)(Point, Point, Point, Point, double);
	double (*func_dA2y_dA1y)(Point, Point, Point, Point, double);
	double (*func_dA2y_dA2x)(Point, Point, Point, Point, double);
	double (*func_dA2y_dA2y)(Point, Point, Point, Point, double);
	double (*func_dA2y_dB1x)(Point, Point, Point, Point, double);
	double (*func_dA2y_dB1y)(Point, Point, Point, Point, double);
	double (*func_dA2y_dB2x)(Point, Point, Point, Point, double);
	double (*func_dA2y_dB2y)(Point, Point, Point, Point, double);
	//
	//	B1x_
	double (*func_dB1x_dA1x)(Point, Point, Point, Point, double);
	double (*func_dB1x_dA1y)(Point, Point, Point, Point, double);
	double (*func_dB1x_dA2x)(Point, Point, Point, Point, double);
	double (*func_dB1x_dA2y)(Point, Point, Point, Point, double);
	double (*func_dB1x_dB1x)(Point, Point, Point, Point, double);
	double (*func_dB1x_dB1y)(Point, Point, Point, Point, double);
	double (*func_dB1x_dB2x)(Point, Point, Point, Point, double);
	double (*func_dB1x_dB2y)(Point, Point, Point, Point, double);
	// B1y_
	double (*func_dB1y_dA1x)(Point, Point, Point, Point, double);
	double (*func_dB1y_dA1y)(Point, Point, Point, Point, double);
	double (*func_dB1y_dA2x)(Point, Point, Point, Point, double);
	double (*func_dB1y_dA2y)(Point, Point, Point, Point, double);
	double (*func_dB1y_dB1x)(Point, Point, Point, Point, double);
	double (*func_dB1y_dB1y)(Point, Point, Point, Point, double);
	double (*func_dB1y_dB2x)(Point, Point, Point, Point, double);
	double (*func_dB1y_dB2y)(Point, Point, Point, Point, double);
	//
	//	B2x_
	double (*func_dB2x_dA1x)(Point, Point, Point, Point, double);
	double (*func_dB2x_dA1y)(Point, Point, Point, Point, double);
	double (*func_dB2x_dA2x)(Point, Point, Point, Point, double);
	double (*func_dB2x_dA2y)(Point, Point, Point, Point, double);
	double (*func_dB2x_dB1x)(Point, Point, Point, Point, double);
	double (*func_dB2x_dB1y)(Point, Point, Point, Point, double);
	double (*func_dB2x_dB2x)(Point, Point, Point, Point, double);
	double (*func_dB2x_dB2y)(Point, Point, Point, Point, double);
	// B2y_
	double (*func_dB2y_dA1x)(Point, Point, Point, Point, double);
	double (*func_dB2y_dA1y)(Point, Point, Point, Point, double);
	double (*func_dB2y_dA2x)(Point, Point, Point, Point, double);
	double (*func_dB2y_dA2y)(Point, Point, Point, Point, double);
	double (*func_dB2y_dB1x)(Point, Point, Point, Point, double);
	double (*func_dB2y_dB1y)(Point, Point, Point, Point, double);
	double (*func_dB2y_dB2x)(Point, Point, Point, Point, double);
	double (*func_dB2y_dB2y)(Point, Point, Point, Point, double);
	
public:
	Constraint();
	~Constraint();

	double	Function();
	double	Derivative(Point *point, bool xy);
	double	SecondDerivative(Point* point_d1, bool xy_d1, Point *point_d2, bool xy_d2);

	void	Fill_Free();
};
