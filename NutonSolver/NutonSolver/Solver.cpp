#define  _CRT_SECURE_NO_WARNINGS

#include "Header.h"

#include <emscripten.h>
#include <emscripten/bind.h>
#include <string>
#include <cstring>

using namespace std;

// расчет системы уравнений
static	double	ft_Calculate_matrix_elem(vector<Constraint> Constraints, vector<Point*> points, double *Ls, size_t equationNmb, size_t derivativeNmb)
{
	double res = 0;

	if (equationNmb < Constraints.size())
	{
		if (derivativeNmb > points.size() * 2)
			return (0.0f);
		return (Constraints[equationNmb].Derivative(points[(derivativeNmb - 1) / 2], derivativeNmb % 2));
	}
	
	if (derivativeNmb - 1 == equationNmb - Constraints.size())
		res += 1;
	
	if (points.size() * 2 >= derivativeNmb)
	{
		for (size_t i = 0; i < Constraints.size(); i++)
		{
			res += Ls[i] * Constraints[i].SecondDerivative(points[(equationNmb - Constraints.size()) / 2], (equationNmb - Constraints.size() + 1) % 2,
				points[(derivativeNmb - 1) / 2], derivativeNmb % 2);
		}
	}
	else
	{
		res += Constraints[derivativeNmb - points.size() * 2 - 1].Derivative(points[(equationNmb - Constraints.size()) / 2],
																				(equationNmb - Constraints.size() + 1) % 2);
	}
	
	return (res);
}

// расчет дельты к координатам
static	double	ft_Calculate_delta(double *lJacobian, vector<Constraint> Constraints, vector<Point*> points, double *Ls)
{
	double	res = 0,
			tmp;

	size_t i = 0;

	for (; i < Constraints.size(); i++)
	{
		res -= Constraints[i].Function() * lJacobian[i];
	}

	for (; i < points.size() * 2 + Constraints.size(); i++)
	{
		if ((i - Constraints.size() + 1) % 2 == 1)
			tmp = points[(i - Constraints.size()) / 2]->dx;
		else
			tmp = points[(i - Constraints.size()) / 2]->dy;

		for (size_t j = 0; j < Constraints.size(); j++)
		{
			tmp += Ls[j] * Constraints[j].Derivative(points[(i - Constraints.size()) / 2], (i - Constraints.size() + 1) % 2);
		}

		res -= tmp * lJacobian[i];
	}

	return (res);
}

static	double	ft_Equation_i(vector<Constraint> Constraints, vector<Point*> points, double* Ls, size_t equationNmb)
{
	double	res = 0;

	if (equationNmb < Constraints.size())
		return (Constraints[equationNmb].Function());

	if ((equationNmb - Constraints.size() + 1) % 2 == 1)
		res = points[(equationNmb - Constraints.size()) / 2]->dx;
	else
		res = points[(equationNmb - Constraints.size()) / 2]->dy;

	for (size_t j = 0; j < Constraints.size(); j++)
	{
		res += Ls[j] * Constraints[j].Derivative(points[(equationNmb - Constraints.size()) / 2], (equationNmb - Constraints.size() + 1) % 2);
	}

	return (res);
}

static	bool	ft_GaussMethod(double** coefficients, double* freeCoefficients, size_t dimension, double* result)
{
	double d, s;

	for (size_t k = 0; k < dimension; k++) // прямой ход
	{
		for (size_t j = k + 1; j < dimension; j++)
		{
			double	a = coefficients[j][k],
					b = coefficients[k][k];
			if (coefficients[k][k] == 0)
				if (MakeDiagonal_NonZero(coefficients, freeCoefficients, 0, dimension) == false)
				{
					coefficients[k][k] = 1;
					//return (false);
				}
			d = coefficients[j][k] / coefficients[k][k]; // формула (1)
			for (size_t i = k; i < dimension; i++)
			{
				coefficients[j][i] = coefficients[j][i] - d * coefficients[k][i]; // формула (2)
			}
			freeCoefficients[j] = freeCoefficients[j] - d * freeCoefficients[k]; // формула (3)
		}
	}

	for (long int k = (long int)dimension - 1; k >= 0; k--) // обратный ход
	{
		d = 0;
		for (long int j = k + 1; j < dimension; j++)
		{
			s = coefficients[k][j] * result[j]; // формула (4)
			d = d + s; // формула (4)
		}
		if (coefficients[k][k] == 0)
			if (MakeDiagonal_NonZero(coefficients, freeCoefficients, 0, dimension) == false)
			{
				//coefficients[k][k] = 1;
				return (false);
			}
		double rrr = coefficients[k][k];
		result[k] = (freeCoefficients[k] - d) / coefficients[k][k]; // формула (4)
	}
	return (true);
}

static	bool	ft_FindDeltas(vector<Constraint> Constraints, vector<Point*> points, double** mJacobian, size_t matrixSize, double* Ls)
{
	double* f;
	double* results;

	if (!(f = (double*)calloc(matrixSize, sizeof(double))))
		return (false);
	if (!(results = (double*)calloc(matrixSize, sizeof(double))))
	{
		free(f);
		return (false);
	}

	// заполнение вектора значений функций на предыдущей итеррации
	for (size_t i = 0; i < matrixSize; i++)
	{
		f[i] = ft_Equation_i(Constraints, points, Ls, i);
		results[i] = 0;
	}

	MakeDiagonal_NonZero(mJacobian, f, 0, matrixSize);

	if (ft_GaussMethod(mJacobian, f, matrixSize, results) == false)
		return (false);

	for (size_t i = 0; (i < points.size() * 2) && (i < matrixSize); i++)
	{
		if ((i + 1) % 2 == 1)
			points[i / 2]->dx -= results[i];
		else
			points[i / 2]->dy -= results[i];
	}
	//cout << endl;
	for (size_t i = points.size() * 2; i < matrixSize; i++)
		Ls[i - points.size() * 2] -= results[i];

	free(results);
	free(f);
	return (true);
}

static	bool	ft_NewtonMethod(vector<Constraint> Constraints, vector<Point*> points)
{
	size_t			matrixSize;
	int				maxIterCount = 75;

	double	* Ls,
			* amendments,
			** mJacobian,
			norm,
			norm_old;

	matrixSize = (Constraints.size() + points.size() * 2);

	Ls = (double*)calloc(points.size() * 2, sizeof(double));
	amendments = (double*)calloc(matrixSize, sizeof(double));
	if (!Ls)			return (false);
	if (!amendments)	return (false);
	for (size_t i = 0; i < points.size() * 2; i++)
		Ls[i] = 1;

	mJacobian = (double**)calloc(matrixSize, sizeof(double*));
	if (!mJacobian)
		return (false);
	for (size_t i = 0; i < matrixSize; i++)
	{
		mJacobian[i] = (double*)calloc(matrixSize, sizeof(double));
		if (!mJacobian[i])
			return (false);
		for (size_t j = 0; j < matrixSize; j++)
			mJacobian[i][j] = 0;
	}

	//-- Calculate first (base) Norm
	norm = 0;
	for (size_t i = 0; i < matrixSize; i++)
	{
		amendments[i] = ft_Equation_i(Constraints, points, Ls, i);
		norm += amendments[i] * amendments[i];
	}
	norm = sqrt(norm);


	do
	{
		norm_old = norm;

		//-- Calculation of the matrix
		for (size_t i = 0; i < matrixSize; i++)
		{
			for (size_t j = 0; j < matrixSize; j++)
			{
				mJacobian[i][j] = ft_Calculate_matrix_elem(Constraints, points, Ls, i, j + 1);
			}
		}

		if (ft_FindDeltas(Constraints, points, mJacobian, matrixSize, Ls) == false)
			return (false);		

		//-- Calculate Norm
		norm = 0;
		for (size_t i = 0; i < matrixSize; i++)
		{
			amendments[i] = ft_Equation_i(Constraints, points, Ls, i);
			norm += amendments[i] * amendments[i];
		}
		norm = sqrt(norm);
		
	} while ((fabs(norm_old - norm) >= eps) && (maxIterCount-- > 0));

	//-- Recalculation of point's coordinates
	for (size_t i = 0; i < points.size(); i++)
	{
		points[i]->x += points[i]->dx;
		points[i]->y += points[i]->dy;

		points[i]->dx = 0;
		points[i]->dy = 0;
	}

	free(Ls);
	free(amendments);
	for (size_t i = 0; i < matrixSize; i++)
		free(mJacobian[i]);
	free(mJacobian);
	return (true);
}

static	void	ft_CycleApproximation(vector<Constraint> Constraints, vector<Point*> points, Point* pointChangeable1, Point* pointChangeable2)
{
	int		maxIterCount = 20;
	Point	pointIdeal1(pointChangeable1->x, pointChangeable1->y),
			pointIdeal2;

	if (pointChangeable2)
		pointIdeal2 = Point(pointChangeable2->x, pointChangeable2->y);

	double norm, norm_old;

	if (ft_NewtonMethod(Constraints, points) == false)
		return;

	norm = pow(pointChangeable1->x - pointIdeal1.x, 2) + pow(pointChangeable1->y - pointIdeal1.y, 2);
	if (pointChangeable2)
		norm += pow(pointChangeable2->x - pointIdeal2.x, 2) + pow(pointChangeable2->y - pointIdeal2.y, 2);
	norm = sqrt(norm);

	do
	{
		norm_old = norm;

		pointChangeable1->x = pointIdeal1.x;
		pointChangeable1->y = pointIdeal1.y;
		if (pointChangeable2)
		{
			pointChangeable2->x = pointIdeal2.x;
			pointChangeable2->y = pointIdeal2.y;
		}

		if (ft_NewtonMethod(Constraints, points) == false)
			return;

		norm = pow(pointChangeable1->x - pointIdeal1.x, 2) + pow(pointChangeable1->y - pointIdeal1.y, 2);
		if (pointChangeable2)
			norm += pow(pointChangeable2->x - pointIdeal2.x, 2) + pow(pointChangeable2->y - pointIdeal2.y, 2);
		norm = sqrt(norm);

	} while ((fabs(norm - norm_old) > eps) && (maxIterCount-- > 0));

}

//extern "C"{
vector<vector<double>>	Solver(string json_str)
{
	Point				*pointChangeable1 = nullptr;
	Point				*pointChangeable2 = nullptr;
	vector<Constraint>	vConstr;
	vector<Point*>		vPointsPtr;
	vector<Point*>		pointsChangeable;

	Json_Read(json_str, &vConstr, &vPointsPtr, &pointChangeable1, &pointChangeable2);

	for (size_t i = 0; i < vPointsPtr.size(); i++)
	{
		if (vPointsPtr[i]->fixed == false)
			pointsChangeable.push_back(vPointsPtr[i]);
	}

	ft_CycleApproximation(vConstr, pointsChangeable, pointChangeable1, pointChangeable2);

	//Json_Write(vPointsPtr);
	vector<double> 			tmpVector;
	vector<vector<double>> 	outputData;		//	Every element is - ID, x, y;

	for (size_t i = 0; i < pointsChangeable.size(); i++)
	{
		tmpVector.push_back((double)pointsChangeable[i]->ID);
		tmpVector.push_back(pointsChangeable[i]->x);
		tmpVector.push_back(pointsChangeable[i]->y);

		outputData.push_back(tmpVector);
		tmpVector.clear();
	}


	for (size_t i = 0; i < vPointsPtr.size(); i++)
		free(vPointsPtr[i]);

	return (outputData);
}
//}

EMSCRIPTEN_BINDINGS(stl_wrappers) {
	emscripten::register_vector<double>("VectorDouble");
	emscripten::register_vector<std::vector<double>>("VectorVectorDouble");
}

EMSCRIPTEN_BINDINGS(module) {
	emscripten::function("Solver", &Solver);
}



//
//int main()
//{
//	Solver(string("{\"Points\":[{\"x\":403.5,\"y\":156.50245154530512,\"id\":2,\"fixed\":false},{\"x\":494.5,\"y\":156.50245154530512,\"id\":3,\"fixed\":false},{\"x\":886.5,\"y\":156.5024515453051,\"id\":4,\"fixed\":false},{\"x\":31.5,\"y\":436.5,\"id\":1,\"fixed\":false}],\"Constraints\":[{\"point1\":2,\"point2\":3,\"point3\":0,\"point4\":0,\"Type\":\"Horizontal_line\",\"value\":0},{\"point1\":3,\"point2\":4,\"point3\":0,\"point4\":0,\"Type\":\"Horizontal_line\",\"value\":0},{\"point1\":3,\"point2\":4,\"point3\":1,\"point4\":2,\"Type\":\"Perpendicularity_of_2_lines\",\"value\":0}],\"MovablePoints_id\":[1,0]}"));
//
//	//Solver(string("{\"Points\":[{\"x\":186.5,\"y\":148.69101063626525,\"id\":1,\"fixed\":false},{\"x\":644.1078155936365,\"y\":148.69101063626525,\"id\":2,\"fixed\":false},{\"x\":675.4993986095989,\"y\":333.50010215221425,\"id\":3,\"fixed\":false},{\"x\":158.44742109835744,\"y\":421.3259427986679,\"id\":4,\"fixed\":false},{\"x\":761.5,\"y\":669.5,\"id\":5,\"fixed\":false},{\"x\":987.5,\"y\":210.5,\"id\":6,\"fixed\":false}],\"Constraints\":[{\"point1\":1,\"point2\":2,\"point3\":0,\"point4\":0,\"Type\":\"Horizontal_line\",\"value\":0},{\"point1\":2,\"point2\":3,\"point3\":3,\"point4\":4,\"Type\":\"Perpendicularity_of_2_lines\",\"value\":0},{\"point1\":3,\"point2\":4,\"point3\":5,\"point4\":6,\"Type\":\"Parallelism_of_2_lines\",\"value\":0}],\"MovablePoints_id\":[5,0]}"));
//
//	//Solver(string("{\"Points\": [{\"x\":239.96297512089856, \"y\" : 368.31119251883246, \"id\" : 2, \"fixed\" : false}, { \"x\":346.0370248791015,\"y\" : 296.6888074811675,\"id\" : 3,\"fixed\" : false }, { \"x\":394.5000971461387,\"y\" : 607.5001438751914,\"id\" : 4,\"fixed\" : false }, { \"x\":544.4393780978064,\"y\" : 506.25946401927314,\"id\" : 5,\"fixed\" : false }] , \"Constraints\" : [{\"point1\":2, \"point2\" : 3, \"point3\" : 4, \"point4\" : 5, \"Type\" : \"Parallelism_of_2_lines\", \"value\" : 0}, { \"point1\":2,\"point2\" : 3,\"point3\" : 0,\"point4\" : 0,\"Type\" : \"Vertical_line\",\"value\" : 0 }] , \"MovablePoints_id\" : [2, 0] }"));
//
//
////	/*Point A1{ 2, 2 };
////	Point A2{ 4, 8 };
////	Point B1{ 3, 2 };
////	Point B2{ 8, 2 };
////
////	vector<Point*> vPointsPtr = { &A1, &A2, &B1, &B2 };
////
////	auto constr = CreateConstraint_Parallelism_of_2_lines(&A1, &A2, &B1, &B2);
////	auto constr_2 = CreateConstraint_Distance_between_2_points(&A1, &A2, 4);
////	auto constr_3 = CreateConstraint_Distance_between_2_points(&B1, &B2, 5);
////	vector<Constraint> vConstr = { constr, constr_2, constr_3 };
////
////
////	ft_CycleApproximation(vConstr, vPointsPtr, &A2, nullptr);
////
////	
////
////	cout << "PointA1 x = " << A1.x << endl;
////	cout << "PointA1 y = " << A1.y << endl << endl;
////	cout << "PointA2 x = " << A2.x << endl;
////	cout << "PointA2 y = " << A2.y << endl << endl;
////	cout << endl;
////	cout << "PointB1 x = " << B1.x << endl;
////	cout << "PointB1 y = " << B1.y << endl << endl;
////	cout << "PointB2 x = " << B2.x << endl;
////	cout << "PointB2 y = " << B2.y << endl << endl;*/
//	return (0);
//}
//
