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
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	x1_h = A1.x + A1.dx + h - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double derivation = (helperFunction(x1_h, x2, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	x1_h = A1.x + A1.dx + h - A2.x - A2.dx - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double derivation = (helperFunction(x1_h, x2, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y1_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1_h = A1.y + A1.dy + h - A2.y - A2.dy;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double derivation = (helperFunction(x1, x2, y1_h, y2) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y1_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1_h = A1.y + A1.dy - A2.y - A2.dy - h;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double derivation = (helperFunction(x1, x2, y1_h, y2) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x2_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	x2_h = B1.x + B1.dx + h - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double derivation = (helperFunction(x1, x2_h, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x2_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	x2_h = B1.x + B1.dx - B2.x - B2.dx - h;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	double derivation = (helperFunction(x1, x2_h, y1, y2) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y2_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	// y1_h = A1.y + A1.dy + h - A2.y - A2.dy;;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;
	y2_h = B1.y + B1.dy + h - B2.y - B2.dy;

	double derivation = (helperFunction(x1, x2, y1, y2_h) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

double Constraint_Angle_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double y2_h;
	double angle;
	// step for counting derivation
	double h = 0.01;

	x1 = A1.x + A1.dx - A2.x - A2.dx;
	// y1_h = A1.y + A1.dy + h - A2.y - A2.dy;;
	y1 = A1.y + A1.dy - A2.y - A2.dy;

	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;
	y2_h = B1.y + B1.dy - B2.y - B2.dy - h;

	double derivation = (helperFunction(x1, x2, y1, y2_h) - helperFunction(x1, x2, y1, y2)) / h;
	return derivation;
}

// set of 2nd derivation for angle constraint
// double Constraint_Angle_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	x1_h = A1.x + A1.dx + h - A2.x - A2.dx;
	x2_h = B1.x + B1.dx + h - B2.x - B2.dx;

	double second_derivation = (helperFunction(x1_h, x2_h, y1, y2) - helperFunction(x1, x2_h, y1, y2) - helperFunction(x1_h, x2, y1, y2)
	+ helperFunction(x1, x2, y1, y2)) / pow(h, 2);

	return second_derivation;
}
double Constraint_Angle_of_2_lines_dA1x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	x1_h = A1.x + A1.dx + h - A2.x - A2.dx;
	x2_h = B1.x + B1.dx - B2.x - B2.dx - h;

	double second_derivation = (helperFunction(x1_h, x2_h, y1, y2) - helperFunction(x1, x2_h, y1, y2) - helperFunction(x1_h, x2, y1, y2)
	+ helperFunction(x1, x2, y1, y2)) / pow(h, 2);

	return second_derivation;
}
// double Constraint_Angle_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	y1_h = A1.y + A1.dy + h - A2.y - A2.dy;
	y2_h = B1.y + B1.dy + h - B2.y - B2.dy;

	double second_derivation = (helperFunction(x1, x2, y1_h, y2_h) - helperFunction(x1,x2,y1, y2_h) - helperFunction(x1,x2,y1_h, y2)
	+ helperFunction(x1,x2,y1,y2)) / pow(h,2);

	return second_derivation;
}

double Constraint_Angle_of_2_lines_dA1y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	y1_h = A1.y + A1.dy + h - A2.y - A2.dy;
	y2_h = B1.y + B1.dy - B2.y - B2.dy - h;

	double second_derivation = (helperFunction(x1, x2, y1_h, y2_h) - helperFunction(x1,x2,y1, y2_h) - helperFunction(x1,x2,y1_h, y2)
	+ helperFunction(x1,x2,y1,y2)) / pow(h,2);

	return second_derivation;
}

// double Constraint_Angle_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	x1_h = A1.x + A1.dx - A2.x - A2.dx - h;
	x2_h = B1.x + B1.dx + h - B2.x - B2.dx;

	double second_derivation = (helperFunction(x1_h, x2_h, y1, y2) - helperFunction(x1, x2_h, y1, y2) - helperFunction(x1_h, x2, y1, y2)
	+ helperFunction(x1, x2, y1, y2)) / pow(h, 2);

	return second_derivation;
}

double Constraint_Angle_of_2_lines_dA2x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	x1_h = A1.x + A1.dx - A2.x - A2.dx - h;
	x2_h = B1.x + B1.dx - B2.x - B2.dx - h;

	double second_derivation = (helperFunction(x1_h, x2_h, y1, y2) - helperFunction(x1, x2_h, y1, y2) - helperFunction(x1_h, x2, y1, y2)
	+ helperFunction(x1, x2, y1, y2)) / pow(h, 2);

	return second_derivation;
}

// double Constraint_Angle_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	y1_h = A1.y + A1.dy - A2.y - A2.dy - h;
	y2_h = B1.y + B1.dy + h - B2.y - B2.dy;

	double second_derivation = (helperFunction(x1, x2, y1_h, y2_h) - helperFunction(x1,x2,y1, y2_h) - helperFunction(x1,x2,y1_h, y2)
	+ helperFunction(x1,x2,y1,y2)) / pow(h,2);

	return second_derivation;
}
double Constraint_Angle_of_2_lines_dA2y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	double x1, x2, y1, y2;
	double x1_h, x2_h , y1_h, y2_h;
	double h = 0.01;
	x1 = A1.x + A1.dx - A2.x - A2.dx;
	y1 = A1.y + A1.dy - A2.y - A2.dy;
	x2 = B1.x + B1.dx - B2.x - B2.dx;
	y2 = B1.y + B1.dy - B2.y - B2.dy;

	y1_h = A1.y + A1.dy - A2.y - A2.dy - h;
	y2_h = B1.y + B1.dy - B2.y - B2.dy - h;

	double second_derivation = (helperFunction(x1, x2, y1_h, y2_h) - helperFunction(x1,x2,y1, y2_h) - helperFunction(x1,x2,y1_h, y2)
	+ helperFunction(x1,x2,y1,y2)) / pow(h,2);

	return second_derivation;
}
// double Constraint_Angle_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle);

// ATTENTION!!!!!!!!!!! Не уверен, что должно быть так, но аналитически выражения для вторых производных получаются одинаковые 
// внезависимости от порядка дифференцирования
double Constraint_Angle_of_2_lines_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dB1x(A1, A2, B1, B2, _angle);
}
double Constraint_Angle_of_2_lines_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2x_dB1x(A1, A2, B1, B2, _angle);
}
// double Constraint_Angle_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dB1y(A1, A2, B1, B2, _angle);
}
double Constraint_Angle_of_2_lines_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2y_dB1y(A1, A2, B1, B2, _angle);
}
// double Constraint_Angle_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1x_dB2x(A1, A2, B1, B2, _angle);
}
double Constraint_Angle_of_2_lines_dB2x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2x_dB2x(A1, A2, B1, B2, _angle);
}
// double Constraint_Angle_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA1y_dB2y(A1, A2, B1, B2, _angle);
}
double Constraint_Angle_of_2_lines_dB2y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle)
{
	return Constraint_Angle_of_2_lines_dA2y_dB2y(A1, A2, B1, B2, _angle);
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

	x1 = B1.x + B1.dx - A1.x - A1.dx;
	y1 = B1.y + B1.dy - A1.y - A1.dy;

	x2 = A2.x + A2.dx - B1.x - B1.dx;
	y2 = A2.y + A2.dy - B1.y - B1.dy;

	return (x1 * x2 + y1 * y2);
}

double Constraint_Belonging_point_to_line_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A2.dx - A2.x + B1.dx + B1.x);
}

double Constraint_Belonging_point_to_line_dA1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A2.dy - A2.y + B1.dy + B1.y);
}

double Constraint_Belonging_point_to_line_dA1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A1.dx - A1.x + B1.dx + B1.x);
}

double Constraint_Belonging_point_to_line_dA2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-A1.dy - A1.y + B1.dy + B1.y);
}

double Constraint_Belonging_point_to_line_dA2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-1);
}

double Constraint_Belonging_point_to_line_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A1.dx + A1.x + A2.dx + A2.x - 2 * B1.dx - 2 * B1.x);
}

double Constraint_Belonging_point_to_line_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dB1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-2);
}

double Constraint_Belonging_point_to_line_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (A1.dy + A1.y + A2.dy + A2.y - 2 * B1.dy - 2 * B1.y);
}

double Constraint_Belonging_point_to_line_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (1);
}

double Constraint_Belonging_point_to_line_dB1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (-2);
}