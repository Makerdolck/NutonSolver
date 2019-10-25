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


// TODO: Добавить выражение 6.45 для вычисления значения вектора на следующей итеррации
static double* ft_JacobiMethod(double** _mJacobian, double* _F, const size_t _n, vector<Point*>& points, vector<Constraint>& vConstr, double* Ls)
{
	// size_t n = sizeof()
	double* y = new double[_n];
	double* f = new double[_n];
	std::vector<double> vDiff;

	// добавить метод для приведения матрицы к ступенчатому виду (возможно)

	// заполнение вектора значений функций на предыдущей итеррации
	for (size_t i = 0; i < _n; i++)
		f[i] = ft_Equation_i(vConstr, points, Ls, i);
	
	// итератор для максимальнго элемента
	auto it_max = vDiff.end();

	do
	{
		for (auto& ptr_it : points)
		{
			ptr_it->x += ptr_it->dx;
			ptr_it->y += ptr_it->dy;
		}
		for (size_t i = 0; i < _n; i++)
		{
			double sum = 0;
			for(size_t j = 0; j < _n; j++)
			{
				if (j != i && (i + 1) % 2 == 1)
					sum += _mJacobian[i][j] * points[i]->dx;
				if (j != i && (i + 1) % 2 == 0)
					sum += _mJacobian[i][j] * points[i]->dy;
			}
			y[i] = (f[i] - sum) / _mJacobian[i][i];

			for (size_t i = 0; i < _n; i++)
				if ((i + 1) % 2)
					vDiff.push_back(y[i] - points[i]->dx);
				else
					vDiff.push_back(y[i] - points[i]->dy);

			// вычисление нормы разности
			it_max = std::max_element(vDiff.begin(), vDiff.end());
		}
	}
	while (*it_max >= eps);
	return y;
}

static void ft_Solver(vector<Constraint> Constraints, vector<Point*> points)
{
	size_t			matrixSize;

	double			*Ls,
					*amendments,
					**mJacobian,
					norm,
					norm_old;

	matrixSize = (Constraints.size() + points.size() * 2);

	Ls = (double*)malloc(points.size() * 2 * sizeof(double));
	amendments = (double*)malloc(matrixSize * sizeof(double));
	if (!Ls)			return;
	if (!amendments)	return;
	for (size_t i = 0; i < points.size() * 2; i++)
		Ls[i] = 1;

	mJacobian = (double**)malloc(matrixSize * sizeof(double*));
	if (!mJacobian)
		return;
	for (size_t i = 0; i < matrixSize; i++)
	{
		mJacobian[i] = (double*)malloc(matrixSize * sizeof(double));
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

		// Inverting_the_matrix(mJacobian, matrixSize); // fix this shit

		//-- Calculation of new deltas
		for (size_t i = 0; (i < points.size() * 2) && (i < matrixSize); i++)
		{
			if ((i + 1) % 2 == 1)
				points[i / 2]->dx += ft_Calculate_delta(mJacobian[i], Constraints, points, Ls);
			else
				points[i / 2]->dy += ft_Calculate_delta(mJacobian[i], Constraints, points, Ls);
		}

		for (size_t i = points.size() * 2; i < Constraints.size(); i++)
			Ls[i - points.size() * 2] += ft_Calculate_delta(mJacobian[i], Constraints, points, Ls);

		//-- Calculate Norm
		norm = 0;
		for (size_t i = 0; i < matrixSize; i++)
		{
			amendments[i] = ft_Equation_i(Constraints, points, Ls, i);
			norm += amendments[i] * amendments[i];
		}
		norm = sqrt(norm);

	} while (fabs(norm_old - norm) >= eps);

	double temp_dx = 0, temp_dy = 0;
	bool isFixed = false;
	//-- Recalculation of point's coordinates
	for (size_t i = 0; i < points.size(); i++)
	{
		// if (points.at(i)->fixed)
		// {
		// 	isFixed = true;
		// 	temp_dx = points.at(i)->dx;
		// 	temp_dy = points.at(i)->dy;
		// }
		// if (!points.at(i)->fixed)
		// {
			points[i]->x += points[i]->dx;
			points[i]->y += points[i]->dy;
		// }
		points[i]->dx = 0;
		points[i]->dy = 0;
	}
	// if (isFixed)
	// {
	// 	auto it = find_if(points.begin(), points.end(), [](Point* a) {return !(a->fixed);});
	// 	(*it)->x -= temp_dx;
	// 	(*it)->y -= temp_dy;
	// }
}

void Solver(vector<Constraint> Constraints, vector<Point*> points, Point *pointChangeable)
{
	Point pointIdeal(pointChangeable->x, pointChangeable->y);

	double norm, norm_old;

	ft_Solver(Constraints, points);
	norm = sqrt(pow(pointChangeable->x - pointIdeal.x, 2) + pow(pointChangeable->y - pointIdeal.y, 2));

	do
	{
		norm_old = norm;

		pointChangeable->x = pointIdeal.x;
		pointChangeable->y = pointIdeal.y;

		ft_Solver(Constraints, points);

		norm = sqrt(pow(pointChangeable->x - pointIdeal.x, 2) + pow(pointChangeable->y - pointIdeal.y, 2));

	} while (fabs(norm - norm_old) > eps);

}

int main()
{
	Point A1{2, 2};
	Point A2{2, 8};
	Point B1{3, 2};
	Point B2{6, 2};

	vector<Point*> vPointsPtr = {&A1, &A2, &B1, &B2};

	auto constr = CreateConstraint_Parallelism_of_2_lines(&A1, &A2, &B1, &B2);
	auto constr_2 = CreateConstraint_Distance_between_2_points(&A1, &A2, 4);
	auto constr_3 = CreateConstraint_Distance_between_2_points(&B1, &B2, 7);
	vector<Constraint> vConstr = {constr, constr_2, constr_3};

	Solver(vConstr, vPointsPtr, &A1);
}