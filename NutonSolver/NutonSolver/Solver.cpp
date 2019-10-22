#define  _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <math.h>
#include <vector>

#include "Header.h"

using namespace std;

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

static	double	ft_Ñalculate_delta(double *lJacobian, vector<Constraint> Constraints, vector<Point*> points, double *Ls)
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

void Solver(vector<Constraint> Constraints, vector<Point*> points)
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

		Inverting_the_matrix(mJacobian, matrixSize);

		//-- Calculation of new deltas
		for (size_t i = 0; (i < points.size() * 2) && (i < matrixSize); i++)
		{
			if ((i + 1) % 2 == 1)
				points[i / 2]->dx += ft_Ñalculate_delta(mJacobian[i], Constraints, points, Ls);
			else
				points[i / 2]->dy += ft_Ñalculate_delta(mJacobian[i], Constraints, points, Ls);
		}

		for (size_t i = points.size() * 2; i < Constraints.size(); i++)
			Ls[i - points.size() * 2] += ft_Ñalculate_delta(mJacobian[i], Constraints, points, Ls);

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
}
