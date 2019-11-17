#include "Matrix.h"
#include <math.h>
#include <stdio.h>

//	
void	Print_matrix(double** arr, size_t size)
{
	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
			printf("%.2f\t", arr[i][j]);
		printf("\n");
	}
}

//	////	///	////

// приведение матрицы коэффициентов к виду с ненулевой диагональю и соответствующее изменение вектора правых частей
// возвращает true - если приведение - успешно
bool	MakeDiagonal_NonZero(double** coefficients, double* rightPart, size_t currColumn, size_t numberOfEquation)
{
	bool result = false;
	size_t i, row = currColumn;
	double tempItem;

	// для матрицы 1x1 ответом является ненулевость ее единственного элемента
	if (currColumn == numberOfEquation - 1) {
		result = coefficients[currColumn][currColumn] != 0;
	}
	else {
		// пока не найдена перестановка приводящая матрицу к виду с ненулевой диагональю и пока мы можем еще просматривать строки ищем перестановку
		while (!result && row < numberOfEquation) {
			// если текущий элемент на диагонали нулевой ищем в столбце дальше ненулевой
			if (coefficients[row][currColumn] != 0) {
				// если элемент ненулевой и не лежит на диагонали меняем местами сотвествующие строки в матрице и элементы в векторе прваых частей
				if (row != currColumn) {
					for (i = 0; i < numberOfEquation; i++) {
						tempItem = coefficients[currColumn][i];
						coefficients[currColumn][i] = coefficients[row][i];
						coefficients[row][i] = tempItem;
					}
					tempItem = rightPart[currColumn];
					rightPart[currColumn] = rightPart[row];
					rightPart[row] = tempItem;
				}
				// рекурсивный вызов фактически для подматрицы с размерностью на 1 меньше
				result = MakeDiagonal_NonZero(coefficients, rightPart, currColumn + 1, numberOfEquation);
				if (result) {
					break;
				}
			}
			row++;
		}
	}

	return result;
}

// было ли найдено решение, если да - итог в параметре solution
int		JacobiMethod(double** coefficients, double* rightPart, size_t numberOfEquation, double* solution, double precision)
{
	bool result;
	int i, j, step = 1;
	double* tempSolution;

	tempSolution = new double[numberOfEquation];

	// приведение матрицы коэффициентов к виду с ненулевой диагональю и соответствующее изменение вектора правых частей
	result = MakeDiagonal_NonZero(coefficients, rightPart, 0, numberOfEquation);

	// если приведение успешно - работаем дальше
	if (result) {
		double fault = precision + 1;

		// преобразуем матрицу и правую часть для дальнейшего решения
		for (i = 0; i < numberOfEquation; i++) {
			for (j = 0; j < numberOfEquation; j++) {
				if (i != j) {
					coefficients[i][j] = -coefficients[i][j] / coefficients[i][i];
				}
			}
			rightPart[i] = rightPart[i] / coefficients[i][i];
			coefficients[i][i] = 0;
		}

		// первое приближение решения - преобразованный вектор правых частей
		for (i = 0; i < numberOfEquation; i++) {
			solution[i] = rightPart[i];
		}

		// пока не найдено решение с допустимй погрешнстью или пока не исчерпан лимит шагов... если все расходится например
		while (fault > precision&& step <= 1000) {

			// поиск новой итерации с ее "самоиспользованием" при расчетах          
			for (j = 0; j < numberOfEquation; j++) {
				tempSolution[j] = 0;
			}
			for (i = 0; i < numberOfEquation; i++) {
				for (j = 0; j < numberOfEquation; j++) {
					tempSolution[i] = tempSolution[i] + coefficients[i][j] * solution[j];
				}
				tempSolution[i] = tempSolution[i] + rightPart[i];
			}

			// расчет погрешности
			fault = 0.0;
			for (j = 0; j < numberOfEquation; j++) {
				fault = fault + (solution[j] - tempSolution[j]) * (solution[j] - tempSolution[j]);
			}
			fault = sqrt(fault);

			// сохранение полученной новой итерации
			for (j = 0; j < numberOfEquation; j++) {
				solution[j] = tempSolution[j];
			}
			step++;
		}
	}
	else {
		step = 1001;
	}


	return step;
}

//	////	///	////

// Matrix transpose function
template <typename T>
static void TransponMtx(T** matr, T** tMatr, int n) {//
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

void Inverting_the_matrix(double** matrix, int size)
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
