#define  _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

#include "Header.h"

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

static	void	ft_Jacobi(size_t matrixSize, double** mJacobian, double* f, double* results)
{
	double* TempX = new double[matrixSize];
	double norm; // норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.

	do {
		for (int i = 0; i < matrixSize; i++) {
			TempX[i] = f[i];
			for (int g = 0; g < matrixSize; g++) {
				if (i != g)
					TempX[i] -= mJacobian[i][g] * results[g];
			}
			TempX[i] /= mJacobian[i][i];
		}
		norm = fabs(results[0] - TempX[0]);
		for (int h = 0; h < matrixSize; h++) {
			if (fabs(results[h] - TempX[h]) > norm)
				norm = fabs(results[h] - TempX[h]);
			results[h] = TempX[h];
		}
	} while (norm > eps);
	delete[] TempX;
}

static	void	ft_GaussMethod(double** coefficients, double* freeCoefficients, size_t dimension, double* result)
{
	double d, s;

	for (size_t k = 0; k < dimension; k++) // прямой ход
	{
		for (size_t j = k + 1; j < dimension; j++)
		{
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
		result[k] = (freeCoefficients[k] - d) / coefficients[k][k]; // формула (4)
	}
}

static	void	ft_FindDeltas(vector<Constraint> Constraints, vector<Point*> points, double** mJacobian, size_t matrixSize, double* Ls)
{
	double* f;
	double* results;

	if (!(f = (double*)calloc(matrixSize, sizeof(double))))
		return;
	if (!(results = (double*)calloc(matrixSize, sizeof(double))))
	{
		free(f);
		return;
	}

	// заполнение вектора значений функций на предыдущей итеррации
	for (size_t i = 0; i < matrixSize; i++)
	{
		f[i] = ft_Equation_i(Constraints, points, Ls, i);
		results[i] = 0;
	}

	MakeDiagonal_NonZero(mJacobian, f, 0, matrixSize);

	ft_GaussMethod(mJacobian, f, matrixSize, results);


	for (size_t i = 0; (i < points.size() * 2) && (i < matrixSize); i++)
	{
		if ((i + 1) % 2 == 1)
			points[i / 2]->dx -= results[i];
		else
			points[i / 2]->dy -= results[i];
	}

	for (size_t i = points.size() * 2; i < matrixSize; i++)
		Ls[i - points.size() * 2] -= results[i];

	free(results);
	free(f);
}

static	void	ft_Solver(vector<Constraint> Constraints, vector<Point*> points)
{
	size_t			matrixSize;

	double* Ls,
		* amendments,
		** mJacobian,
		norm,
		norm_old;

	matrixSize = (Constraints.size() + points.size() * 2);

	Ls = (double*)calloc(points.size() * 2, sizeof(double));
	amendments = (double*)calloc(matrixSize, sizeof(double));
	if (!Ls)			return;
	if (!amendments)	return;
	for (size_t i = 0; i < points.size() * 2; i++)
		Ls[i] = 1;

	mJacobian = (double**)calloc(matrixSize, sizeof(double*));
	if (!mJacobian)
		return;
	for (size_t i = 0; i < matrixSize; i++)
	{
		mJacobian[i] = (double*)calloc(matrixSize, sizeof(double));
		if (!mJacobian[i])
			return;
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

		ft_FindDeltas(Constraints, points, mJacobian, matrixSize, Ls);

		//-- Calculate Norm
		norm = 0;
		for (size_t i = 0; i < matrixSize; i++)
		{
			amendments[i] = ft_Equation_i(Constraints, points, Ls, i);
			norm += amendments[i] * amendments[i];
		}
		norm = sqrt(norm);

	} while (fabs(norm_old - norm) >= eps);

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
}

void Solver(vector<Constraint> Constraints, vector<Point*> points, Point* pointChangeable1, Point* pointChangeable2)
{
	Point	pointIdeal1(pointChangeable1->x, pointChangeable1->y),
		pointIdeal2;

	if (pointChangeable2)
		pointIdeal2 = Point(pointChangeable2->x, pointChangeable2->y);

	double norm, norm_old;

	ft_Solver(Constraints, points);

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

		ft_Solver(Constraints, points);

		norm = pow(pointChangeable1->x - pointIdeal1.x, 2) + pow(pointChangeable1->y - pointIdeal1.y, 2);
		if (pointChangeable2)
			norm += pow(pointChangeable2->x - pointIdeal2.x, 2) + pow(pointChangeable2->y - pointIdeal2.y, 2);
		norm = sqrt(norm);

	} while (fabs(norm - norm_old) > eps);

}

int main()
{
	Point A1{ 2, 2 };
	Point A2{ 4, 8 };
	Point B1{ 3, 2 };
	Point B2{ 8, 2 };

	vector<Point*> vPointsPtr = { &A1, &A2, &B1, &B2 };

	auto constr = CreateConstraint_Parallelism_of_2_lines(&A1, &A2, &B1, &B2);
	auto constr_2 = CreateConstraint_Distance_between_2_points(&A1, &A2, 4);
	auto constr_3 = CreateConstraint_Distance_between_2_points(&B1, &B2, 5);
	vector<Constraint> vConstr = { constr, constr_2, constr_3 };

	Solver(vConstr, vPointsPtr, &A2, nullptr);

	cout << "PointA1 x = " << A1.x << endl;
	cout << "PointA1 y = " << A1.y << endl << endl;
	cout << "PointA2 x = " << A2.x << endl;
	cout << "PointA2 y = " << A2.y << endl << endl;
	cout << endl;
	cout << "PointB1 x = " << B1.x << endl;
	cout << "PointB1 y = " << B1.y << endl << endl;
	cout << "PointB2 x = " << B2.x << endl;
	cout << "PointB2 y = " << B2.y << endl << endl;
}