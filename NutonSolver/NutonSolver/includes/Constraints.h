#pragma once

#include "Point.h"

/*																							Constraint_Match_2_points_x				*/

double Constraint_Match_2_points_x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Match_2_points_x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Match_2_points_x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);

/*																							Constraint_Match_2_points_y				*/

double Constraint_Match_2_points_y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Match_2_points_y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Match_2_points_y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);

/*																							Constraint_Distance_between_2_points	*/

double Constraint_Distance_between_2_points(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dA1x(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dA1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Distance_between_2_points_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Distance_between_2_points_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Distance_between_2_points_dA1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Distance_between_2_points_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Distance_between_2_points_dB1x(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dB1x_dB1x(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dB1y(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double distance);
double Constraint_Distance_between_2_points_dB1y_dB1y(Point A1, Point A2, Point B1, Point B2, double distance);

/*																							Constraint_Parallelism_of_2_lines		*/

double Constraint_Parallelism_of_2_lines(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA1x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA1x_dB2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA1y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA1y_dB2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA2x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA2x_dB2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA2y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dA2y_dB2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB1x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB1y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB2x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB2x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB2y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Parallelism_of_2_lines_dB2y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);

/*																							Constraint_Perpendicularity_of_2_lines	*/

double Constraint_Perpendicularity_of_2_lines(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA1x_dB2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA1y_dB2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA2x_dB2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dA2y_dB2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB2x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Perpendicularity_of_2_lines_dB2y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);

/*																							Constraint_Angle_of_2_lines				*/

double Constraint_Angle_of_2_lines(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA1y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2x_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dA2y_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2x_dA2x(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _angle);
double Constraint_Angle_of_2_lines_dB2y_dA2y(Point A1, Point A2, Point B1, Point B2, double _angle);

/*																							Constraint_Horizontal_line				*/

double Constraint_Horizontal_line(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Horizontal_line_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Horizontal_line_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);

/*																							Constraint_Vertical_line				*/

double Constraint_Vertical_line(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Vertical_line_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Vertical_line_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);

/*																							Constraint_Belonging_point_to_line		*/

double Constraint_Belonging_point_to_line(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1x_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1x_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1x_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1y_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1y_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1y_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);

//	---	---	--- New
double Constraint_Belonging_point_to_line_dA1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);

double Constraint_Belonging_point_to_line_dA1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA1y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);

double Constraint_Belonging_point_to_line_dA2x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2x_dB1y(Point A1, Point A2, Point B1, Point B2, double _value);

double Constraint_Belonging_point_to_line_dA2y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dA2y_dB1x(Point A1, Point A2, Point B1, Point B2, double _value);

double Constraint_Belonging_point_to_line_dB1x_dA1y(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1x_dA2y(Point A1, Point A2, Point B1, Point B2, double _value);

double Constraint_Belonging_point_to_line_dB1y_dA1x(Point A1, Point A2, Point B1, Point B2, double _value);
double Constraint_Belonging_point_to_line_dB1y_dA2x(Point A1, Point A2, Point B1, Point B2, double _value);