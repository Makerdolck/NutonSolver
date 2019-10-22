#pragma once

#include "Point.h"
#include "Constraint.h"
#include "Constraints.h"

Constraint	CreateConstraint_Match_2_points_x(Point* A1, Point* B1);
Constraint	CreateConstraint_Match_2_points_y(Point* A1, Point* B1);
Constraint	CreateConstraint_Distance_between_2_points(Point* A1, Point* B1, double distance);
Constraint	CreateConstraint_Parallelism_of_2_lines(Point* A1, Point* A2, Point* B1, Point* B2);
Constraint	CreateConstraint_Perpendicularity_of_2_lines(Point* A1, Point* A2, Point* B1, Point* B2);
Constraint	CreateConstraint_Horizontal_line(Point* A1, Point* A2);
Constraint	CreateConstraint_Vertical_line(Point* A1, Point* A2);
Constraint	CreateConstraint_Belonging_point_to_line(Point* A1, Point* A2, Point* B1);


