#include <cmath>
#include "Constraints.h"


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

//double Constraint_Angle_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle)
//{
//	double x1, x2, y1, y2;
//	double res, f, g, df, dg;
//
//	x1 = A1.x + A1.dx - A2.x - A2.dx;
//	y1 = A1.y + A1.dy - A2.y - A2.dy;
//
//	x2 = B1.x + B1.dx - B2.x - B2.dx;
//	y2 = B1.y + B1.dy - B2.y - B2.dy;
//
//	f = (x1 * x2 + y1 * y2);
//	g = sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2);
//
//	df = B1.q + B1.x - B2.q - B2.x;
//	dg =							;	// (f�g)' = f'�g + f�g'
//
//	res = df * g - f * dg / (pow(g, 2));
//
//	return (res);
//}

/*																							Constraint_Horizontal_line				*/

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
