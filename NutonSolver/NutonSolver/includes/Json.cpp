#include <iomanip>
#include <iostream>
#include <string>
#include "Json.h"

static	Point		*ft_GetPointByID(vector<Point*> *points, size_t id)
{
	for (size_t i = 0; i < points->size(); i++)
	{
		if (id == (*points)[i]->ID)
			return ((*points)[i]);
	}
	return (nullptr);
}

static	Constraint	ft_Determine_Constraint(string type, vector<Point*>* points_All, size_t	points_id[4], double value)
{
	Constraint		constraint;
	vector<Point*>	points;

	for (int i = 0; i < 4; i++)
	{
		if (points_id[i] == 0)
			continue;
		points.push_back(ft_GetPointByID(points_All, points_id[i]));
	}

	if (type == string("Distance_between_2_points"))
		return (CreateConstraint_Distance_between_2_points(points[0], points[1], value));
	if (type == string("Parallelism_of_2_lines"))
		return (CreateConstraint_Parallelism_of_2_lines(points[0], points[1], points[2], points[3]));
	if (type == string("Perpendicularity_of_2_lines"))
		return (CreateConstraint_Perpendicularity_of_2_lines(points[0], points[1], points[2], points[3]));
	if (type == string("Horizontal_line"))
		return (CreateConstraint_Horizontal_line(points[0], points[1]));
	if (type == string("Vertical_line"))
		return (CreateConstraint_Vertical_line(points[0], points[1]));
	if (type == string("Belonging_point_to_line"))
		return (CreateConstraint_Belonging_point_to_line(points[0], points[1], points[2]));
	if (type == string("Angle_between_2_lines"))
		return (CreateConstraint_Angle_between_2_lines(points[0], points[1], points[2], points[3], value));
	exit(1);
}

bool				Json_Read(	string				json_str,
								vector<Constraint>	*Constraints, 
								vector<Point*>		*points, 
								Point				**pointChangeable1, 
								Point				**pointChangeable2)
{
	size_t			points_id[4];
	string			constraint_type;
	json_space		json,
					js_group;
	Point			*point;
	Constraint		constraint;
	std::ofstream	outputFile;
	std::ifstream	inputFile;

	outputFile.open("_TMP_.json", ios_base::trunc);
	if (!outputFile.is_open())
		return (false);
	outputFile << json_str;
	outputFile.close();

	inputFile.open("_TMP_.json");
	if (!inputFile.is_open())
		return (false);

	inputFile >> json;
	inputFile.close();

	js_group = json["Points"];
	for (long int i = 0; i < js_group.end() - js_group.begin(); i++)
	{
		point = new Point();
		point->x		= js_group[i]["x"];
		point->y		= js_group[i]["y"];
		point->ID		= js_group[i]["id"];
		point->fixed	= js_group[i]["fixed"];

		points->push_back(point);
	}

	js_group = json["Constraints"];
	for (long int i = 0; i < js_group.end() - js_group.begin(); i++)
	{
		constraint_type = js_group[i]["Type"];

		points_id[0] = js_group[i]["point1"];
		points_id[1] = js_group[i]["point2"];
		points_id[2] = js_group[i]["point3"];
		points_id[3] = js_group[i]["point4"];
		
		Constraints->push_back(ft_Determine_Constraint(	constraint_type,
														points,
														points_id,
														js_group[i]["value"]));
	}

	*pointChangeable1 = ft_GetPointByID(points, json["MovablePoints_id"][0]);
	
	if (json["MovablePoints_id"][1] != 0)
		*pointChangeable2 = ft_GetPointByID(points, json["MovablePoints_id"][1]);

	return (true);
}


bool				Json_Write(vector<Point*> points)
{
	json_space		json, 
					point_js;

	std::ofstream	file_out("calculated.json");
	if (!file_out.is_open())
		return (false);

	for (size_t i = 0; i < points.size(); i++)
	{
		point_js["id"] = points[i]->ID;
		point_js["x"] = points[i]->x;
		point_js["y"] = points[i]->y;

		json.push_back(point_js);
	}

	file_out << std::setw(4) << json << endl;

	file_out.close();

	return (true);
}