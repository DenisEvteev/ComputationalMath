#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <cassert>
#include <fstream>
#include <cstring>

using namespace std;

void find_L(double **L, double **A, int size);
double sumSquareElemStr(int line, double **L);
double sumLineRowMultip(int line, double **L);
void readMatrixFromFile(ifstream& in, int size, double **A);
void printMatrix(int size, double **A);
double **transpose(int size, double **M);
void remove(int size, double **M);
double **init(int size);
double *createFreeRow(int size, double **A);
double *additional_vector(int size, double *f, double **L);
double sumStrVec(double **L, double *v, int low, int up);
double *solution(int size, double *v, double **L_t);
double norm_difference(int size, double* v);

int main(int argc, char **argv)
{
	if (argc != 2) {
		cerr << "Give me a filename to read the matrix, please!" << endl;
		exit(EXIT_FAILURE);
	}
	int size = 0; //will represent an exponent of the initial matrix we will read it from a file
	ifstream stream(argv[1]);

	if (!stream.is_open()) {
		cerr << "Cannot open file " << argv[1] << endl;
		exit(EXIT_FAILURE);
	}
	stream >> size;
	assert(size > 0); // read exponent of a matrix

	try {
		double **A = init(size);
		readMatrixFromFile(stream, size, A);
		double *b = createFreeRow(size, A);
		
		double **L = init(size);

		find_L(L, A, size);

		double **L_t = transpose(size, L);
		
		double *v = additional_vector(size, b, L);

		double *u = solution(size, v, L_t);
		// cout << "Result:\n";
		// for (int i = 0; i < size; ++i)
		// 	cout << u[i] << endl;

		// cout << "The norm of difference of X and Y :\n";
		cout << norm_difference(size, u) << endl; // we'd like to calculate the norm
		
		//clear resources
		stream.close();
		remove(size, A);
		delete [] b;
		delete [] v;
		delete [] u;
		remove(size, L);
		remove(size, L_t);
	
	}catch(const exception& ex) {
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	return 0;
}

double norm_difference(int size, double *u)
{
	double *tmp = new double[size];
	memset(tmp, 0.0, size);
	
	for (int i = 0; i < size; ++i) {
		tmp[i] = (i + 1) - u[i];
	}
	//find norm of the vector tmp
	double square_sum = 0;
	for (int i = 0; i < size; ++i) {
		square_sum += tmp[i] * tmp[i];
	}

	delete [] tmp;
	return sqrt(square_sum);
}

double *solution(int size, double *v, double **L_t)
{
	double *u = new double[size];
	memset(u, 0.0, size);
	for (int i = size - 1; i >= 0; --i) {
		if (L_t[i][i] == 0.0)
			throw runtime_error("Division by zero!");
		u[i] = (v[i] - sumStrVec(L_t, u, i + 1, size)) / L_t[i][i];
	}
	return u;
}

double sumStrVec(double **L, double *v, int low, int up)
{
	double sum = 0;
	for(int ip_ = low; ip_ < up; ++ip_) {
		
		sum += L[(low == 0) ? up : (low - 1)][ip_] * v[ip_];
	}
	return sum;
}

double *additional_vector(int size, double *f, double **L)
{
	double *v = new double[size];
	memset(v, 0.0, size);
	for (int i = 0; i < size; ++i) {
		if (L[i][i] == 0.0) // TODO: create more sophisticated check
			throw runtime_error("Division by zero!");
		v[i] = (f[i] - sumStrVec(L, v, 0, i)) / L[i][i];
	}
	return v;
}

double *createFreeRow(int size, double **A)
{
	double *f = new double[size];
	memset(f, 0.0, size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			f[i] += A[i][j] * (j + 1);
		}
	}
	return f;
}

double **init(int size)
{
	double **M = new double*[size];
	for (int i = 0; i < size; ++i) {
		M[i] = new double[size];
		memset(M[i], 0.0, size);
	}
	return M;
}

void remove(int size, double **M)
{
	for (int i = 0; i < size; ++i)
		delete [] M[i];
	delete[] M;
}

void printMatrix(int size, double **A)
{
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

void readMatrixFromFile(ifstream& in, int size, double **A)
{
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			in >> A[i][j];
		}
	}
}

void find_L(double **L, double **A, int size)
{
	double diag = 0;
	for (int str = 0; str < size; ++str) {
		diag = A[str][str] - sumSquareElemStr(str, L);
		if (diag < 0)
			throw runtime_error("Negative value under the root!");
		L[str][str] = sqrt(diag);

		if (L[str][str] == 0)
			throw runtime_error("Division by zero!");

		for (int down = str + 1; down < size; ++down) {
			L[down][str] = (A[str][down] - sumLineRowMultip(str, L)) / L[str][str];
		}
	}
}

double sumLineRowMultip(int line, double **L)
{
	double sum = 0;
	for (int ip_ = 0; ip_ < line; ++ip_) {
		sum += (L[line][ip_] * L[line + 1][ip_]);
	}
	return sum;
}

double sumSquareElemStr(int line, double **L)
{
	double sum = 0;
	for (int ip_ = 0; ip_ < line; ++ip_) {
		sum += L[line][ip_] * L[line][ip_];
	}
	return sum;
}

double **transpose(int size, double **M)
{
	double **M_t = init(size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			M_t[i][j] = M[j][i];
		}

	}
	return M_t;
}
