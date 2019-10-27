
#include "CreateConstraint.h"
#include "Point.h"

#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "nlohmann/json.hpp"

using json_space = nlohmann::json;
using namespace std;

bool	Json_Read(	string				json_str, 
					vector<Constraint>	*Constraints, 
					vector<Point*>		*points, 
					Point				**pointChangeable1, 
					Point				**pointChangeable2);

bool	Json_Write(	vector<Point*>		points);