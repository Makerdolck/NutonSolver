#include "CreateConstraint.h"

Constraint	CreateConstraint_Match_2_points_x(Point *A1, Point *B1)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointB1 = B1;

	constr.func = &Constraint_Match_2_points_x;
	constr.func_dA1x = &Constraint_Match_2_points_x_dA1x;
	constr.func_dB1x = &Constraint_Match_2_points_x_dB1x;

	return (constr);
}

Constraint	CreateConstraint_Match_2_points_y(Point* A1, Point* B1)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointB1 = B1;

	constr.func = &Constraint_Match_2_points_y;
	constr.func_dA1y = &Constraint_Match_2_points_y_dA1y;
	constr.func_dB1y = &Constraint_Match_2_points_y_dB1y;

	return (constr);
}

Constraint	CreateConstraint_Distance_between_2_points(Point* A1, Point* B1, double distance)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointB1 = B1;
	constr.value   = distance;

	constr.func = &Constraint_Distance_between_2_points;

	constr.func_dA1x = &Constraint_Distance_between_2_points_dA1x;
	constr.func_dA1x_dA1x = &Constraint_Distance_between_2_points_dA1x_dA1x;
	constr.func_dA1x_dB1x = &Constraint_Distance_between_2_points_dA1x_dB1x;

	constr.func_dA1y = &Constraint_Distance_between_2_points_dA1y;
	constr.func_dA1y_dA1y = &Constraint_Distance_between_2_points_dA1y_dA1y;
	constr.func_dA1y_dB1y = &Constraint_Distance_between_2_points_dA1y_dB1y;

	constr.func_dB1x = &Constraint_Distance_between_2_points_dB1x;
	constr.func_dB1x_dA1x = &Constraint_Distance_between_2_points_dB1x_dA1x;
	constr.func_dB1x_dB1x = &Constraint_Distance_between_2_points_dB1x_dB1x;

	constr.func_dB1y = &Constraint_Distance_between_2_points_dB1y;
	constr.func_dB1y_dA1y = &Constraint_Distance_between_2_points_dB1y_dA1y;
	constr.func_dB1y_dB1y = &Constraint_Distance_between_2_points_dB1y_dB1y;

	return (constr);
}

Constraint	CreateConstraint_Parallelism_of_2_lines(Point* A1, Point* A2, Point* B1, Point* B2)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointA2 = A2;
	constr.pointB1 = B1;
	constr.pointB2 = B2;

	constr.func = &Constraint_Parallelism_of_2_lines;
	constr.func_dA1x = &Constraint_Parallelism_of_2_lines_dA1x;
	constr.func_dA1x_dB1y = &Constraint_Parallelism_of_2_lines_dA1x_dB1y;
	constr.func_dA1x_dB2y = &Constraint_Parallelism_of_2_lines_dA1x_dB2y;

	constr.func_dA1y = &Constraint_Parallelism_of_2_lines_dA1y;
	constr.func_dA1y_dB1x = &Constraint_Parallelism_of_2_lines_dA1y_dB1x;
	constr.func_dA1y_dB2x = &Constraint_Parallelism_of_2_lines_dA1y_dB2x;

	constr.func_dA2x = &Constraint_Parallelism_of_2_lines_dA2x;
	constr.func_dA2x_dB1y = &Constraint_Parallelism_of_2_lines_dA2x_dB1y;
	constr.func_dA2x_dB2y = &Constraint_Parallelism_of_2_lines_dA2x_dB2y;

	constr.func_dA2y = &Constraint_Parallelism_of_2_lines_dA2y;
	constr.func_dA2y_dB1x = &Constraint_Parallelism_of_2_lines_dA2y_dB1x;
	constr.func_dA2y_dB2x = &Constraint_Parallelism_of_2_lines_dA2y_dB2x;

	constr.func_dB1x = &Constraint_Parallelism_of_2_lines_dB1x;
	constr.func_dB1x_dA1y = &Constraint_Parallelism_of_2_lines_dB1x_dA1y;
	constr.func_dB1x_dA2y = &Constraint_Parallelism_of_2_lines_dB1x_dA2y;

	constr.func_dB1y = &Constraint_Parallelism_of_2_lines_dB1y;
	constr.func_dB1y_dA1x = &Constraint_Parallelism_of_2_lines_dB1y_dA1x;
	constr.func_dB1y_dA2x = &Constraint_Parallelism_of_2_lines_dB1y_dA2x;

	constr.func_dB2x = &Constraint_Parallelism_of_2_lines_dB2x;
	constr.func_dB2x_dA1y = &Constraint_Parallelism_of_2_lines_dB2x_dA1y;
	constr.func_dB2x_dA2y = &Constraint_Parallelism_of_2_lines_dB2x_dA2y;

	constr.func_dB2y = &Constraint_Parallelism_of_2_lines_dB2y;
	constr.func_dB2y_dA1x = &Constraint_Parallelism_of_2_lines_dB2y_dA1x;
	constr.func_dB2y_dA2x = &Constraint_Parallelism_of_2_lines_dB2y_dA2x;

	return (constr);
}

Constraint	CreateConstraint_Perpendicularity_of_2_lines(Point* A1, Point* A2, Point* B1, Point* B2)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointA2 = A2;
	constr.pointB1 = B1;
	constr.pointB2 = B2;

	constr.func = &Constraint_Perpendicularity_of_2_lines;
	constr.func_dA1x = &Constraint_Perpendicularity_of_2_lines_dA1x;
	constr.func_dA1x_dB1x = &Constraint_Perpendicularity_of_2_lines_dA1x_dB1x;
	constr.func_dA1x_dB2x = &Constraint_Perpendicularity_of_2_lines_dA1x_dB2x;

	constr.func_dA1y = &Constraint_Perpendicularity_of_2_lines_dA1y;
	constr.func_dA1y_dB1y = &Constraint_Perpendicularity_of_2_lines_dA1y_dB1y;
	constr.func_dA1y_dB2y = &Constraint_Perpendicularity_of_2_lines_dA1y_dB2y;

	constr.func_dA2x = &Constraint_Perpendicularity_of_2_lines_dA2x;
	constr.func_dA2x_dB1x = &Constraint_Perpendicularity_of_2_lines_dA2x_dB1x;
	constr.func_dA2x_dB2x = &Constraint_Perpendicularity_of_2_lines_dA2x_dB2x;

	constr.func_dA2y = &Constraint_Perpendicularity_of_2_lines_dA2y;
	constr.func_dA2y_dB1y = &Constraint_Perpendicularity_of_2_lines_dA2y_dB1y;
	constr.func_dA2y_dB2y = &Constraint_Perpendicularity_of_2_lines_dA2y_dB2y;

	constr.func_dB1x = &Constraint_Perpendicularity_of_2_lines_dB1x;
	constr.func_dB1x_dA1x = &Constraint_Perpendicularity_of_2_lines_dB1x_dA1x;
	constr.func_dB1x_dA2x = &Constraint_Perpendicularity_of_2_lines_dB1x_dA2x;

	constr.func_dB1y = &Constraint_Perpendicularity_of_2_lines_dB1y;
	constr.func_dB1y_dA1y = &Constraint_Perpendicularity_of_2_lines_dB1y_dA1y;
	constr.func_dB1y_dA2y = &Constraint_Perpendicularity_of_2_lines_dB1y_dA2y;

	constr.func_dB2x = &Constraint_Perpendicularity_of_2_lines_dB2x;
	constr.func_dB2x_dA1x = &Constraint_Perpendicularity_of_2_lines_dB2x_dA1x;
	constr.func_dB2x_dA2x = &Constraint_Perpendicularity_of_2_lines_dB2x_dA2x;

	constr.func_dB2y = &Constraint_Perpendicularity_of_2_lines_dB2y;
	constr.func_dB2y_dA1y = &Constraint_Perpendicularity_of_2_lines_dB2y_dA1y;
	constr.func_dB2y_dA2y = &Constraint_Perpendicularity_of_2_lines_dB2y_dA2y;

	return (constr);
}

Constraint	CreateConstraint_Horizontal_line(Point* A1, Point* A2)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointA2 = A2;

	constr.func = &Constraint_Horizontal_line;
	constr.func_dA2y = &Constraint_Horizontal_line_dA2y;
	constr.func_dA1y = &Constraint_Horizontal_line_dA1y;

	return (constr);
}

Constraint	CreateConstraint_Vertical_line(Point* A1, Point* A2)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointA2 = A2;

	constr.func = &Constraint_Vertical_line;
	constr.func_dA2x = &Constraint_Vertical_line_dA2x;
	constr.func_dA1x = &Constraint_Vertical_line_dA1x;

	return (constr);
}

Constraint	CreateConstraint_Belonging_point_to_line(Point* A1, Point* A2, Point *B1)
{
	Constraint constr;

	constr.pointA1 = A1;
	constr.pointA2 = A2;
	constr.pointB1 = B1;

	constr.func = &Constraint_Belonging_point_to_line;
	constr.func_dA1x = &Constraint_Belonging_point_to_line_dA1x;
	constr.func_dA1x_dA2x = &Constraint_Belonging_point_to_line_dA1x_dA2x;
	constr.func_dA1x_dB1x = &Constraint_Belonging_point_to_line_dA1x_dB1x;
	constr.func_dA1x_dA2y = &Constraint_Belonging_point_to_line_dA1x_dA2y;
	constr.func_dA1x_dB1y = &Constraint_Belonging_point_to_line_dA1x_dB1y;

	constr.func_dA1y = &Constraint_Belonging_point_to_line_dA1y;
	constr.func_dA1y_dA2y = &Constraint_Belonging_point_to_line_dA1y_dA2y;
	constr.func_dA1y_dB1y = &Constraint_Belonging_point_to_line_dA1y_dB1y;
	constr.func_dA1y_dA2x = &Constraint_Belonging_point_to_line_dA1y_dA2x;
	constr.func_dA1y_dB1x = &Constraint_Belonging_point_to_line_dA1y_dB1x;

	constr.func_dA2x = &Constraint_Belonging_point_to_line_dA2x;
	constr.func_dA2x_dA1x = &Constraint_Belonging_point_to_line_dA2x_dA1x;
	constr.func_dA2x_dB1x = &Constraint_Belonging_point_to_line_dA2x_dB1x;
	constr.func_dA2x_dA1y = &Constraint_Belonging_point_to_line_dA2x_dA1y;
	constr.func_dA2x_dB1y = &Constraint_Belonging_point_to_line_dA2x_dB1y;

	constr.func_dA2y = &Constraint_Belonging_point_to_line_dA2y;
	constr.func_dA2y_dA1y = &Constraint_Belonging_point_to_line_dA2y_dA1y;
	constr.func_dA2y_dB1y = &Constraint_Belonging_point_to_line_dA2y_dB1y;
	constr.func_dA2y_dA1x = &Constraint_Belonging_point_to_line_dA2y_dA1x;
	constr.func_dA2y_dB1x = &Constraint_Belonging_point_to_line_dA2y_dB1x;

	constr.func_dB1x = &Constraint_Belonging_point_to_line_dB1x;
	constr.func_dB1x_dA1x = &Constraint_Belonging_point_to_line_dB1x_dA1x;
	constr.func_dB1x_dA2x = &Constraint_Belonging_point_to_line_dB1x_dA2x;
	constr.func_dB1x_dB1x = &Constraint_Belonging_point_to_line_dB1x_dB1x;
	constr.func_dB1x_dA1y = &Constraint_Belonging_point_to_line_dB1x_dA1y;
	constr.func_dB1x_dA2y = &Constraint_Belonging_point_to_line_dB1x_dA2y;

	constr.func_dB1y = &Constraint_Belonging_point_to_line_dB1y;
	constr.func_dB1y_dA1y = &Constraint_Belonging_point_to_line_dB1y_dA1y;
	constr.func_dB1y_dA2y = &Constraint_Belonging_point_to_line_dB1y_dA2y;
	constr.func_dB1y_dB1y = &Constraint_Belonging_point_to_line_dB1y_dB1y;
	constr.func_dB1y_dA1x = &Constraint_Belonging_point_to_line_dB1y_dA1x;
	constr.func_dB1y_dA2x = &Constraint_Belonging_point_to_line_dB1y_dA2x;

	//// Point projection to Line

	//double	lambda,
	//	x = A2->x - A,
	//	y = Vector.Y;
	//// TODO : Normalization

	//lambda = (-1) * (x * A1->x + y * A1->y  +
	//		x * (-1) * (point->X) + y * (-1) * (point->Y)) /
	//		(x * x + y * y);

	//B1->x = A1->x + x * lambda;
	//B1->y = A1->y + y * lambda;

	return (constr);
}

Constraint CreateConstraint_Angle_between_2_lines(Point* A1, Point* A2, Point* B1, Point* B2, double _angle)
{
	Constraint constr;
	constr.pointA1 = A1;
	constr.pointA2 = A2;
	constr.pointB1 = B1;
	constr.pointB2 = B2;

	constr.func = &Constraint_Angle_of_2_lines;

	constr.func_dA1x = &Constraint_Angle_of_2_lines_dA1x;

	constr.func_dA1x_dA1x = &Constraint_Angle_of_2_lines_dA1x_dA1x;
	constr.func_dA1x_dA1y = &Constraint_Angle_of_2_lines_dA1x_dA1y;
	constr.func_dA1x_dA2x = &Constraint_Angle_of_2_lines_dA1x_dA2x;
	constr.func_dA1x_dA2y = &Constraint_Angle_of_2_lines_dA1x_dA2y;
	constr.func_dA1x_dB1x = &Constraint_Angle_of_2_lines_dA1x_dB1x;
	constr.func_dA1x_dB1y = &Constraint_Angle_of_2_lines_dA1x_dB1y;
	constr.func_dA1x_dB2x = &Constraint_Angle_of_2_lines_dA1x_dB2x;
	constr.func_dA1x_dB2y = &Constraint_Angle_of_2_lines_dA1x_dB2y;

	constr.func_dA1y = &Constraint_Angle_of_2_lines_dA1y;

	constr.func_dA1y_dA1x = &Constraint_Angle_of_2_lines_dA1y_dA1x;
	constr.func_dA1y_dA1y = &Constraint_Angle_of_2_lines_dA1y_dA1y;
	constr.func_dA1y_dA2x = &Constraint_Angle_of_2_lines_dA1y_dA2x;
	constr.func_dA1y_dA2y = &Constraint_Angle_of_2_lines_dA1y_dA2y;
	constr.func_dA1y_dB1x = &Constraint_Angle_of_2_lines_dA1y_dB1x;
	constr.func_dA1y_dB1y = &Constraint_Angle_of_2_lines_dA1y_dB1y;
	constr.func_dA1y_dB2x = &Constraint_Angle_of_2_lines_dA1y_dB2x;
	constr.func_dA1y_dB2y = &Constraint_Angle_of_2_lines_dA1y_dB2y;

	constr.func_dA2x = &Constraint_Angle_of_2_lines_dA2x;

	constr.func_dA2x_dA1x = &Constraint_Angle_of_2_lines_dA2x_dA1x;
	constr.func_dA2x_dA1y = &Constraint_Angle_of_2_lines_dA2x_dA1y;
	constr.func_dA2x_dA2x = &Constraint_Angle_of_2_lines_dA2x_dA2x;
	constr.func_dA2x_dA2y = &Constraint_Angle_of_2_lines_dA2x_dA2y;
	constr.func_dA2x_dB1x = &Constraint_Angle_of_2_lines_dA2x_dB1x;
	constr.func_dA2x_dB1y = &Constraint_Angle_of_2_lines_dA2x_dB1y;
	constr.func_dA2x_dB2x = &Constraint_Angle_of_2_lines_dA2x_dB2x;
	constr.func_dA2x_dB2y = &Constraint_Angle_of_2_lines_dA2x_dB2y;

	constr.func_dA2y = &Constraint_Angle_of_2_lines_dA2y;

	constr.func_dA2y_dA1x = &Constraint_Angle_of_2_lines_dA2y_dA1x;
	constr.func_dA2y_dA1y = &Constraint_Angle_of_2_lines_dA2y_dA1y;
	constr.func_dA2y_dA2x = &Constraint_Angle_of_2_lines_dA2y_dA2x;
	constr.func_dA2y_dA2y = &Constraint_Angle_of_2_lines_dA2y_dA2y;
	constr.func_dA2y_dB1x = &Constraint_Angle_of_2_lines_dA2y_dB1x;
	constr.func_dA2y_dB1y = &Constraint_Angle_of_2_lines_dA2y_dB1y;
	constr.func_dA2y_dB2x = &Constraint_Angle_of_2_lines_dA2y_dB2x;
	constr.func_dA2y_dB2y = &Constraint_Angle_of_2_lines_dA2y_dB2y;

	constr.func_dB1x = &Constraint_Angle_of_2_lines_dB1x;

	constr.func_dB1x_dA1x = &Constraint_Angle_of_2_lines_dB1x_dA1x;
	constr.func_dB1x_dA1y = &Constraint_Angle_of_2_lines_dB1x_dA1y;
	constr.func_dB1x_dA2x = &Constraint_Angle_of_2_lines_dB1x_dA2x;
	constr.func_dB1x_dA2y = &Constraint_Angle_of_2_lines_dB1x_dA2y;
	constr.func_dB1x_dB1x = &Constraint_Angle_of_2_lines_dB1x_dB1x;
	constr.func_dB1x_dB1y = &Constraint_Angle_of_2_lines_dB1x_dB1y;
	constr.func_dB1x_dB2x = &Constraint_Angle_of_2_lines_dB1x_dB2x;
	constr.func_dB1x_dB2y = &Constraint_Angle_of_2_lines_dB1x_dB2y;

	constr.func_dB1y = &Constraint_Angle_of_2_lines_dB1y;

	constr.func_dB1y_dA1x = &Constraint_Angle_of_2_lines_dB1y_dA1x;
	constr.func_dB1y_dA1y = &Constraint_Angle_of_2_lines_dB1y_dA1y;
	constr.func_dB1y_dA2x = &Constraint_Angle_of_2_lines_dB1y_dA2x;
	constr.func_dB1y_dA2y = &Constraint_Angle_of_2_lines_dB1y_dA2y;
	constr.func_dB1y_dB1x = &Constraint_Angle_of_2_lines_dB1y_dB1x;
	constr.func_dB1y_dB1y = &Constraint_Angle_of_2_lines_dB1y_dB1y;
	constr.func_dB1y_dB2x = &Constraint_Angle_of_2_lines_dB1y_dB2x;
	constr.func_dB1y_dB2y = &Constraint_Angle_of_2_lines_dB1y_dB2y;

	constr.func_dB2x = &Constraint_Angle_of_2_lines_dB2x;

	constr.func_dB2x_dA1x = &Constraint_Angle_of_2_lines_dB2x_dA1x;
	constr.func_dB2x_dA1y = &Constraint_Angle_of_2_lines_dB2x_dA1y;
	constr.func_dB2x_dA2x = &Constraint_Angle_of_2_lines_dB2x_dA2x;
	constr.func_dB2x_dA2y = &Constraint_Angle_of_2_lines_dB2x_dA2y;
	constr.func_dB2x_dB1x = &Constraint_Angle_of_2_lines_dB2x_dB1x;
	constr.func_dB2x_dB1y = &Constraint_Angle_of_2_lines_dB2x_dB1y;
	constr.func_dB2x_dB2x = &Constraint_Angle_of_2_lines_dB2x_dB2x;
	constr.func_dB2x_dB2y = &Constraint_Angle_of_2_lines_dB2x_dB2y;

	constr.func_dB2y = &Constraint_Angle_of_2_lines_dB2y;

	constr.func_dB2y_dA1x = &Constraint_Angle_of_2_lines_dB2y_dA1x;
	constr.func_dB2y_dA1y = &Constraint_Angle_of_2_lines_dB2y_dA1y;
	constr.func_dB2y_dA2x = &Constraint_Angle_of_2_lines_dB2y_dA2x;
	constr.func_dB2y_dA2y = &Constraint_Angle_of_2_lines_dB2y_dA2y;
	constr.func_dB2y_dB1x = &Constraint_Angle_of_2_lines_dB2y_dB1x;
	constr.func_dB2y_dB1y = &Constraint_Angle_of_2_lines_dB2y_dB1y;
	constr.func_dB2y_dB2x = &Constraint_Angle_of_2_lines_dB2y_dB2x;
	constr.func_dB2y_dB2y = &Constraint_Angle_of_2_lines_dB2y_dB2y;

	

	return constr;
}
