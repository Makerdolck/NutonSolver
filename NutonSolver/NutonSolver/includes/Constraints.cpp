#include <cmath>
#include "Constraints.h"

// helping functions for angle constraint between two lines
double helperFunction(double x1, double x2, double y1, double y2)
{

	return (x1 * x2 + y1 * y2) / ((sqrt(x1 * x1 + y1 * y1) * (sqrt(x2 * x2 + y2 * y2))));
}


/*																							Constraint_Match_2_points_x				*/

double Constraint_Match_2_points_x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return A1.dx - B1.dx;
}

double Constraint_Match_2_points_x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return 1;
}

double Constraint_Match_2_points_x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return -1;
}

/*																							Constraint_Match_2_points_y				*/

double Constraint_Match_2_points_y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return A1.dy - B1.dy;
}

double Constraint_Match_2_points_y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return 1;
}

double Constraint_Match_2_points_y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return -1;
}

/*																							Constraint_Distance_between_2_points	*/

double Constraint_Distance_between_2_points(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (pow((B1.x + B1.dx - A1.x - A1.dx), 2) + pow((B1.y + B1.dy - A1.y - A1.dy), 2) - pow(distance, 2));
}

double Constraint_Distance_between_2_points_dA1x(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (-2 * (-A1.x + B1.x - A1.dx + B1.dx));
}

double Constraint_Distance_between_2_points_dA1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (2);
}

double Constraint_Distance_between_2_points_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-2);
}

double Constraint_Distance_between_2_points_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-2 * (-A1.y + B1.y - A1.dy + B1.dy));
}

double Constraint_Distance_between_2_points_dA1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (2);
}

double Constraint_Distance_between_2_points_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-2);
}

double Constraint_Distance_between_2_points_dB1x(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (2 * (-A1.x + B1.x - A1.dx + B1.dx));
}

double Constraint_Distance_between_2_points_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (-2);
}

double Constraint_Distance_between_2_points_dB1x_dB1x(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (2);
}

double Constraint_Distance_between_2_points_dB1y(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (2 * (-A1.y + B1.y - A1.dy + B1.dy));
}

double Constraint_Distance_between_2_points_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (-2);
}

double Constraint_Distance_between_2_points_dB1y_dB1y(Point A1, Point A2, Point B1, Point B2, double distance)
{
	return (2);
}

/*																							Constraint_Parallelism_of_2_lines		*/

double Constraint_Parallelism_of_2_lines(Point A1, Point A2, Point B1, Point B2, double _value)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;


	//// Normalization
	//double length = sqrt(pow(x1, 2) + pow(y1, 2));
	//x1 /= length;
	//y1 /= length;

	///*length = sqrt(pow(x2, 2) + pow(y2, 2));
	//x2 /= length;
	//y2 /= length;*/

	return (x1 * y2 - x2 * y1);
}

double Constraint_Parallelism_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (B1.dy + B1.y - B2.dy - B2.y);
}

double Constraint_Parallelism_of_2_lines_dA1x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dA1x_dB2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-B1.dx - B1.x + B2.dx + B2.x);
}

double Constraint_Parallelism_of_2_lines_dA1y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dA1y_dB2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-B1.dy - B1.y + B2.dy + B2.y);
}

double Constraint_Parallelism_of_2_lines_dA2x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dA2x_dB2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (B1.dx + B1.x - B2.dx - B2.x);
}

double Constraint_Parallelism_of_2_lines_dA2y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dA2y_dB2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A1.dy - A1.y + A2.dy + A2.y);
}

double Constraint_Parallelism_of_2_lines_dB1x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dB1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A1.dx + A1.x - A2.dx - A2.x);
}

double Constraint_Parallelism_of_2_lines_dB1y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dB1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A1.dy + A1.y - A2.dy - A2.y);
}

double Constraint_Parallelism_of_2_lines_dB2x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Parallelism_of_2_lines_dB2x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A1.dx - A1.x + A2.dx + A2.x);
}

double Constraint_Parallelism_of_2_lines_dB2y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Parallelism_of_2_lines_dB2y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

/*																							Constraint_Perpendicularity_of_2_lines	*/

double Constraint_Perpendicularity_of_2_lines(Point A1, Point A2, Point B1, Point B2, double _value)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;


	//// Normalization
	//double length = sqrt(pow(x1, 2) + pow(y1, 2));
	//x1 /= length;
	//y1 /= length;

	///*length = sqrt(pow(x2, 2) + pow(y2, 2));
	//x2 /= length;
	//y2 /= length;*/

	return (x1 * x2 + y1 * y2);
}

double Constraint_Perpendicularity_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (B1.dx + B1.x - B2.dx - B2.x);
}

double Constraint_Perpendicularity_of_2_lines_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dA1x_dB2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (B1.dy + B1.y - B2.dy - B2.y);
}

double Constraint_Perpendicularity_of_2_lines_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dA1y_dB2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-B1.dx - B1.x + B2.dx + B2.x);
}

double Constraint_Perpendicularity_of_2_lines_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dA2x_dB2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-B1.dy - B1.y + B2.dy + B2.y);
}

double Constraint_Perpendicularity_of_2_lines_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dA2y_dB2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A1.dx + A1.x - A2.dx - A2.x);
}

double Constraint_Perpendicularity_of_2_lines_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A1.dy + A1.y - A2.dy - A2.y);
}

double Constraint_Perpendicularity_of_2_lines_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A1.dx - A1.x + A2.dx + A2.x);
}

double Constraint_Perpendicularity_of_2_lines_dB2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dB2x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Perpendicularity_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A1.dy - A1.y + A2.dy + A2.y);
}

double Constraint_Perpendicularity_of_2_lines_dB2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Perpendicularity_of_2_lines_dB2y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

/*																							Constraint_Angle_of_2_lines				*/

double Constraint_Angle_of_2_lines(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double angle;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	angle = (x1 * x2 + y1 * y2) / ((sqrt(x1 * x1 + y1 * y1) * (sqrt(x2 * x2 + y2 * y2))));

	return (angle - _angle);
}

double Constraint_Angle_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	x1_h = A1.x + A1.dx + h - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double first_m = B1.x - B2.x + B1.dx - B2.dx / 
	(sqrt(pow(A1.y - A2.y + A1.dy - A2.dy, 2) + pow(A1.dx + A1.x - A2.x - A2.dx, 2)) * sqrt(pow(B1.x - B2.x + B1.dx - B2.dx, 2) + 
			pow(B1.y - B1.dy - B2.y - B2.dy, 2)));

	double second_m = (A1.dx + A1.x - A2.dx - A2.x) * ((A1.dx + A1.x - A2.dx - A2.x) * (B1.x - B2.x + B1.dx - B2.dx) +
	(A1.y - A2.y + A1.dy - A2.dy) * (B1.y - B1.dy - B2.y - B2.dy));

	double second_mm = pow(pow(A1.y - A2.y + A1.dy - A2.dy, 2) + pow(A1.dx + A1.x - A2.x - A2.dx, 2), 1.5) * 
	sqrt(pow(B1.x - B2.x + B1.dx - B2.dx, 2) + pow(B1.y - B1.dy - B2.y - B2.dy, 2));

	// double derivation = (helperFunction(x1_h, x2, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m / second_mm;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	// double x1_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = x1 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2))); 
	double second_m = x2 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));

	// double derivation = (helperFunction(x1_h, x2, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y1_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = y2 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = y1 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));

	// double derivation = (helperFunction(x1, x2, y1_h, y2) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y1_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = y1 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(y2, 2) + pow(x1, 2)));
	double second_m = y2 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));

	// double derivation = (helperFunction(x1, x2, y1_h, y2) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x2_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = x1 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(y2, 2) + pow(x2, 2)));
	double second_m = x2 * (x1 * x2 + y1 * y2) / (sqrt(pow(y1, 2) + pow(x1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));

	// double derivation = (helperFunction(x1, x2_h, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x2_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = x2 * (x1 * x2 + y1 * y2) / (sqrt(pow(y1, 2) + pow(x1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double second_m = x1 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));

	// double derivation = (helperFunction(x1, x2_h, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y2_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = y1 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = y2 * (x1 * x2 + y1 * y2) / (sqrt(pow(y1, 2) + pow(x1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));

	// double derivation = (helperFunction(x1, x2, y1, y2_h) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y2_h;
	//double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = y2 * (x1 * x2 + y1 * y2) / (sqrt(pow(y1, 2) + pow(x1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double second_m = y1 / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));

	// double derivation = (helperFunction(x1, x2, y1, y2_h) - helperFunction(x1, x2, y1, y2)) / h;
	double derivation = first_m - second_m;
	return derivation;
}

// set of 2nd derivation for angle constraint
// double Constraint_Angle_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = 3 * pow(x1, 2) * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = 2 * x1 * x2 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double third_m = (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));

	double derivation = first_m - second_m - third_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA1x_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = y1 * x2 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = x1 * y2 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double third_m = 3 * x1 * y1 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));

	double derivation = -first_m - second_m + third_m;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = 3 * pow(x1, 2) * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = 2 * x1 * x2 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double third_m = (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));

	return -first_m + second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = (y1 * x2 + x1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = 3 * x1 * y1 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));


	return first_m - second_m;
}

double Constraint_Angle_of_2_lines_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = pow(x1, 2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = x1 * x2 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double third_m = (1 - pow(x2, 2)) / (sqrt(pow(y1, 2) + pow(x1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));

	return -first_m + second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA1x_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = x1 * y1 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = x2 * y2 / (sqrt(pow(y1, 2) + pow(x1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double third_m = x1 * y2 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * pow(pow(x2, 2) + pow(y2, 2), 1.5));


	return -first_m - second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA1x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = pow(x1, 2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = x1 * x2 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double third_m = (1 - pow(x2, 2)) / (sqrt(pow(x1, 2) + pow(y1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));

	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dA1x_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = x1 * y1 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = x2 * y2 / (sqrt(pow(y1, 2) + pow(x1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double third_m = x1 * y2 * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * pow(pow(x2, 2) + pow(y2, 2), 1.5));


	return first_m + second_m - third_m;
}
// double Constraint_Angle_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle);

double Constraint_Angle_of_2_lines_dA1y_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	// double x1, x2, y1, y2;

	// x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// // x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	// y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	// x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	// y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h
	return Constraint_Angle_of_2_lines_dA1x_dA1y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dA1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = 3 * pow(y1, 2) * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = 2 * y1 * x2 / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double third_m = (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));

	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dA1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = (y1 * x2 + x1 * y2) / (pow(pow(x1, 2) + pow(y1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = 3 * x1 * y1 * (x1 * x2 + y1 * y2) / (pow(pow(x1, 2) + pow(y1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));


	return first_m - second_m;
}

double Constraint_Angle_of_2_lines_dA1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = 3 * pow(y1, 2) * (x1 * x2 + y1 * y2) / (pow(pow(y1, 2) + pow(x1, 2), 2.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = (2 * y1 * y2 + (x1 * x2 + y1 * y2)) / (pow(pow(y1, 2) + pow(x1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));


	return -first_m + second_m;
}

double Constraint_Angle_of_2_lines_dA1y_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = x1 * y1 / (pow(pow(x1, 2) + pow(y1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = y1 * x2 * (x1 * x2 + y1 * y2) / (pow(pow(x1, 2) + pow(y1, 2), 1.5) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double third_m = x2 * y2 / (sqrt(pow(x1, 2) + pow(y1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));

	return -first_m + second_m - third_m;
}

double Constraint_Angle_of_2_lines_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double first_m = pow(y1, 2) / (pow(pow(x1, 2) + pow(y1, 2), 1.5) * sqrt(pow(x2, 2) + pow(y2, 2)));
	double second_m = y1 * y2 * (x1 * x2 + y1 * y2) / (pow(pow(x1, 2) + pow(y1, 2), 1.5) * pow(pow(x2,2) + pow(y2, 2), 1.5));
	double third_m = pow(y2, 2) / (sqrt(pow(x1, 2) + pow(y1, 2)) * pow(pow(x2, 2) + pow(y2, 2), 1.5));
	double fourth_m = 1 / (sqrt(pow(x1, 2) + pow(y1, 2)) * sqrt(pow(x2, 2) + pow(y2, 2)));


	return -first_m + second_m - third_m + fourth_m;
}

double Constraint_Angle_of_2_lines_dA1y_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = x1 * y1 / (pow(length1, 1.5) * sqrt(length2));
	double second_m = y1 * x2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = x2 * y2 / (sqrt(length1) * pow(length2, 1.5));

	return first_m - second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA1y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = pow(y1, 2) / (pow(length1, 1.5) * sqrt(length2));
	double second_m = y1 * y2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = pow(y2, 2) / (sqrt(length1) * pow(length2, 2));
	double fourth_m = 1 / (sqrt(length1) * sqrt(length2));

	return first_m - second_m + third_m - fourth_m;
}

// TODO: dA2x, dA2y, dB1x, dB1y, dB2x, dB2y

double Constraint_Angle_of_2_lines_dA2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dA2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dA2x_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dA2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dA2x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * pow(x1, 2) * (x1 * x2 + y1 * y2) / (pow(length1, 2.5) * sqrt(length2));
	double second_m = 2 * x1 * x2 / (pow(length1, 1.5) * sqrt(length2));
	double third_m = (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * sqrt(length2));


	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dA2x_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = y1 * x2 / (pow(length1, 1.5) * sqrt(length2));
	double second_m = x1 * y2 / (pow(length1, 1.5) * sqrt(length2));
	double third_m = 3 * x1 * y1 * (x1 * x2 + y1 * y2) / (pow(length1, 2.5) * sqrt(length2));


	return -first_m - second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = pow(x1, 2) / (pow(length1, 1.5) * sqrt(length2));
	double second_m = x1 * x2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = 1 / (sqrt(length1) * sqrt(length2));
	double fourth_m = pow(x2, 2) / (sqrt(length1) * pow(length2, 1.5));


	return first_m - second_m - third_m + fourth_m;
}

// HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double Constraint_Angle_of_2_lines_dA2x_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = x1 * y1 / (pow(length1, 1.5) * sqrt(length2));
	double second_m = x2 * y2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = x1 * y2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * (pow(length2, 1.5)));


	return first_m + second_m - third_m;
}

double Constraint_Angle_of_2_lines_dA2x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = pow(x1, 2) / (pow(length1, 1.5) * sqrt(length2));
	double second_m = x1 * x2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = 1 / (sqrt(length1) * sqrt(length2));
	double fourth_m = pow(x2, 2) / (sqrt(length1) * pow(length2, 1.5));


	return -first_m + second_m + third_m - fourth_m;
}

double Constraint_Angle_of_2_lines_dA2x_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = x1 * y1 / (pow(length1, 1.5) * sqrt(length2));
	double second_m = x2 * y2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = x1 * y2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));


	return -first_m - second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA2y_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	// double x1, x2, y1, y2;

	// x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// // x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	// y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	// x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	// y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	// double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	// double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	
	return Constraint_Angle_of_2_lines_dA1x_dA2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dA2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	// double x1, x2, y1, y2;

	// x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// // x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	// y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	// x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	// y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	// double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	// double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	
	return Constraint_Angle_of_2_lines_dA1y_dA2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dA2y_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	// double x1, x2, y1, y2;

	// x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// // x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	// y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	// x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	// y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	// double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	// double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	
	return Constraint_Angle_of_2_lines_dA2x_dA2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dA2y_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = x1 * y1 / (pow(length1, 1.5) * sqrt(length2));
	double second_m = y1 * x2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = x2 * y2 / (sqrt(length1) * pow(length2, 1.5));

	return first_m - second_m + third_m;
}

double Constraint_Angle_of_2_lines_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = pow(y1, 2) / (pow(length1, 1.5) * sqrt(length2));
	double second_m = y1 * y2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = pow(y2, 2) / (sqrt(length1) * pow(length2, 1.5));
	double fourth_m = 1 / (sqrt(length1) * sqrt(length2));

	return first_m - second_m + third_m - fourth_m;
}

double Constraint_Angle_of_2_lines_dA2y_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = x1 * y1 / (pow(length1, 1.5) * sqrt(length2));
	double second_m = y1 * x2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = x2 * y2 / (sqrt(length1) * pow(length2, 1.5));


	return -first_m + second_m - third_m;
}

double Constraint_Angle_of_2_lines_dA2y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = pow(y1, 2) / (pow(length1, 1.5) * sqrt(length2));
	double second_m = y1 * y2 * (x1 * x2 + y1 * y2) / (pow(length1, 1.5) * pow(length2, 1.5));
	double third_m = pow(y2, 2) / (sqrt(length1) * pow(length2, 1.5));
	double fourth_m = 1 / (sqrt(length1) * sqrt(length2));


	return -first_m + second_m - third_m + fourth_m;
}

double Constraint_Angle_of_2_lines_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dB1x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1x_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dB1x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2x_dB1x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2y_dB1x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * pow(x2, 2) * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = 2 * x1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 1.5));


	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dB1x_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * x2 * y2 * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = y1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = x1 * y2 / (sqrt(length1) * pow(length2, 1.5));


	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dB1x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * pow(x2, 2) * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = 2 * x1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 1.5));


	return -first_m + second_m + third_m;
}

double Constraint_Angle_of_2_lines_dB1x_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * x2 * y2 * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = y1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = x1 * y2 / (sqrt(length1) * pow(length2, 1.5));


	return -first_m + second_m + third_m;
}

double Constraint_Angle_of_2_lines_dB1y_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dB1y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dB1y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2x_dB1y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2y_dB1y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1y_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dB1x_dB1y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * pow(y2, 2) * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = 2 * y1 * y2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m  = (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 1.5));


	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dB1y_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * x2 * y2 * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = y1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = x1 * y2 / (sqrt(length1) * pow(length2, 1.5));


	return -first_m + second_m + third_m;
}

double Constraint_Angle_of_2_lines_dB1y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка


	double first_m = 3 * pow(y2, 2) * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = 2 * y1 * y2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 1.5));	


	return -first_m + second_m + third_m;
}

double Constraint_Angle_of_2_lines_dB2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dB2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2x_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dB2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2x_dB2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2x_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2y_dB2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dB1x_dB2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2x_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dB1y_dB2x(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * pow(x2, 2) * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = 2 * x1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 1.5));


	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dB2x_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * x2 * y2 * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = y1 * x2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = x1 * y2 / (sqrt(length1) * pow(length2, 1.5));

	return first_m - second_m - third_m;
}

double Constraint_Angle_of_2_lines_dB2y_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2x_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2y_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dB1x_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dB1y_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dB2x_dB2y(A1, A2, B1, B2, _angle);
}

double Constraint_Angle_of_2_lines_dB2y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;

	x1 = A1.x + A1.dx - A2.x - A2.dx;	// a_1 + a - a_3 - c
	// x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;	// a_2 + b - a_4 - d

	x2 = B1.x + B1.dx - B2.x - B2.dx;	// b_1 + e - b_3 - g
	y2 = B1.y + B1.dy - B2.y - B2.dy;	// b_2 + f - b_4 - h

	double length1 = pow(x1, 2) + pow(y1, 2);		// квадрат длины первого отрезка
	double length2 = pow(x2, 2) + pow(y2, 2);		// квадрат длины второго отрезка

	double first_m = 3 * pow(y2, 2) * (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 2.5));
	double second_m = 2 * y1 * y2 / (sqrt(length1) * pow(length2, 1.5));
	double third_m = (x1 * x2 + y1 * y2) / (sqrt(length1) * pow(length2, 1.5));


	return first_m - second_m - third_m;
}

double Constraint_Horizontal_line(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A2.y + A2.dy - A1.y - A1.dy);
}

double Constraint_Horizontal_line_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Horizontal_line_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

/*																							Constraint_Vertical_line				*/

double Constraint_Vertical_line(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A2.x + A2.dx - A1.x - A1.dx);
}

double Constraint_Vertical_line_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Vertical_line_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

/*																							Constraint_Belonging_point_to_line		*/

double Constraint_Belonging_point_to_line(Point A1, Point A2, Point B1, Point B2, double _value)
{
	double x1, x2, y1, y2;

	/*x1 = B1.x + B1.dx - A1.x - A1.dx;
	y1 = B1.y + B1.dy - A1.y - A1.dy;

	x2 = A2.x + A2.dx - B1.x - B1.dx;
	y2 = A2.y + A2.dy - B1.y - B1.dy;*/

	//return (x1 * x2 + y1 * y2);	//x1 * x2 + y1 * y2

	x1 = B1.x + B1.dx - A1.x - A1.dx;
	y1 = B1.y + B1.dy - A1.y - A1.dy;

	x2 = A2.x + A2.dx - A1.x - A1.dx;
	y2 = A2.y + A2.dy - A1.y - A1.dy;

	// Normalization
	double length = sqrt(pow(x1, 2) + pow(y1, 2));
	//x1 /= length;
	//y1 /= length;

	length = sqrt(pow(x2, 2) + pow(y2, 2));
	x2 /= length;
	y2 /= length;

	return (x1 * y2 - x2 * y1);
}

double Constraint_Belonging_point_to_line_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-A2.dx - A2.x + B1.dx + B1.x);
	return (-A2.dy - A2.y + B1.dy + B1.y);
}

double Constraint_Belonging_point_to_line_dA1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-1);
	return (0);
}

double Constraint_Belonging_point_to_line_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

//	---	---	--- New
double Constraint_Belonging_point_to_line_dA1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dA1x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}
//	---	---	---

double Constraint_Belonging_point_to_line_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-A2.dy - A2.y + B1.dy + B1.y);
	return (-B1.dx - B1.x + A2.dx + A2.x);
}

double Constraint_Belonging_point_to_line_dA1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-1);
	return (0);
}

double Constraint_Belonging_point_to_line_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

//	---	---	--- New
double Constraint_Belonging_point_to_line_dA1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dA1y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}
//	---	---	---

double Constraint_Belonging_point_to_line_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-A1.dx - A1.x + B1.dx + B1.x);
	return (A1.dy + A1.y - B1.dy - B1.y);
}

double Constraint_Belonging_point_to_line_dA2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-1);
	return (0);
}

double Constraint_Belonging_point_to_line_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

//	---	---	--- New
double Constraint_Belonging_point_to_line_dA2x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dA2x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}
//	---	---	---

double Constraint_Belonging_point_to_line_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-A1.dy - A1.y + B1.dy + B1.y);
	return (-A1.dx - A1.x + B1.dx + B1.x);
}

double Constraint_Belonging_point_to_line_dA2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-1);
	return (0);
}

double Constraint_Belonging_point_to_line_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

//	---	---	--- New
double Constraint_Belonging_point_to_line_dA2y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dA2y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}
//	---	---	---

double Constraint_Belonging_point_to_line_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (A1.dx + A1.x + A2.dx + A2.x - 2 * B1.dx - 2 * B1.x);
	return (A2.dy + A2.y - A1.dy - A1.y);
}

double Constraint_Belonging_point_to_line_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

double Constraint_Belonging_point_to_line_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

double Constraint_Belonging_point_to_line_dB1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-2);
	return (0);
}

//	---	---	--- New
double Constraint_Belonging_point_to_line_dB1x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dB1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}
//	---	---	---

double Constraint_Belonging_point_to_line_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (A1.dy + A1.y + A2.dy + A2.y - 2 * B1.dy - 2 * B1.y);
	return (A1.dx + A1.x - A2.dx - A2.x);
}

double Constraint_Belonging_point_to_line_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

double Constraint_Belonging_point_to_line_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (1);
	return (0);
}

double Constraint_Belonging_point_to_line_dB1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	//return (-2);
	return (0);
}

//	---	---	--- New
double Constraint_Belonging_point_to_line_dB1y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dB1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}
//	---	---	---