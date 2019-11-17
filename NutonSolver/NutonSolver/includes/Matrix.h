#pragma once

#include <cstring>

int		JacobiMethod(double** coefficients, double* rightPart, size_t numberOfEquation, double* solution, double precision);
bool	MakeDiagonal_NonZero(double** coefficients, double* rightPart, size_t currColumn, size_t numberOfEquation);

void	Inverting_the_matrix(double** matrix, int size);
void	Print_matrix(double** arr, size_t size);