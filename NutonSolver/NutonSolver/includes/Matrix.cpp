#include "Matrix.h"
#include <math.h>
#include <stdio.h>

//	
void	Print_matrix(double** arr, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			printf("%.3f\t", arr[i][j]);
		printf("\n");
	}
}

//	////	///	////

// Matrix transpose function
template <typename T>
static	void	TransponMtx(T** matr, T** tMatr, int n) {//
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			tMatr[j][i] = matr[i][j];
}

// Function to release the memory
template <typename T>
static	void	FreeMem(T** matr, int n)
{
	for (int i = 0; i < n; i++)
		delete[] matr[i];
	delete[] matr;
}

// A function of deleting a row and column
static	void	Get_matr(double** matr, int n, double** temp_matr, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matr[ki][kj] = matr[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}

// The function to compute the determinant of the matrix
static	double	Det(double** matr, int n)
{
	double temp = 0;	// temporary variable to store the determinant
	int k = 1;			// rank
	if (n < 1) {
		printf("Invalid matrix size!!!\n");
		return 0;
	}
	else if (n == 1)
		temp = matr[0][0];
	else if (n == 2)
		temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
	else {
		for (int i = 0; i < n; i++) {
			int m = n - 1;
			double** temp_matr = new double* [m];
			for (int j = 0; j < m; j++)
				temp_matr[j] = new double[m];
			Get_matr(matr, n, temp_matr, 0, i);
			temp = temp + k * matr[0][i] * Det(temp_matr, m);
			k = -k;
			FreeMem(temp_matr, m);
		}
	}
	return temp;
}

//--	Inverting

void			Inverting_the_matrix(double** matrix, int size)
{
	double det;

	double** obr_matr = new double* [size];
	double** tobr_matr = new double* [size];
	for (int i = 0; i < size; i++) {
		obr_matr[i] = new double[size];
		tobr_matr[i] = new double[size];
	}

	det = Det(matrix, size);
	if (det) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				int m = size - 1;
				double** temp_matr = new double* [m];
				for (int k = 0; k < m; k++)
					temp_matr[k] = new double[m];
				Get_matr(matrix, size, temp_matr, i, j);
				obr_matr[i][j] = pow(-1.0, (double)i + j + 2) * Det(temp_matr, m) / det;
				FreeMem(temp_matr, m);
			}
		}
	}
	else
		return;
	// Transposing the matrix
	TransponMtx(obr_matr, tobr_matr, size);
	//FreeMem(tobr_matr, size);
	FreeMem(obr_matr, size);

	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; j++)
			matrix[i][j] = tobr_matr[i][j];
	}
}
