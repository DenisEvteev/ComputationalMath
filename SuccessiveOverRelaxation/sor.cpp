#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <stdexcept>

using matrix_t = std::vector<std::vector<double>>;
using column_t = std::vector<double>;

#define ERROR (1e-4)

const double W = 1.5;

#define PRINT_MATRIX(A)                      \
do{                                          \
	for (int i = 0; i < A.size(); ++i) {     \
		for (int j = 0; j < A.size(); ++j) { \
			std::cout << A[i][j] << " ";     \
		}                                    \
		std::cout << std::endl;              \
	}                                        \
}while(0)

#define PRINT_COLUMN(X)                      \
do{                                          \
	for (int i = 0; i < X.size(); ++i) {     \
		std::cout << X[i] << " ";            \
	}                                        \
	std::cout << std::endl;                  \
}while(0)

void read_matrix_from_file(std::ifstream& in, matrix_t& M);
std::vector<double> find_solution(const matrix_t& A, const std::vector<double>& f);
double euclidean_norm_residual(const matrix_t& A, const column_t& f, const column_t& x);
double previous_sum(const matrix_t& A, const column_t& cur, int row);
double right_part(const matrix_t& A, const column_t& f, const column_t& x, int row);

int main(int argc, char **argv)
{
	if (argc != 2) {
		std::cerr << "bad amount of arguments\n";
		std::exit(EXIT_FAILURE);
	}
	std::ifstream stream(argv[1]);

	if (!stream.is_open()) {
		std::cerr << "Cannot open file " << argv[1] << std::endl;
		std::exit(EXIT_FAILURE);
	}
	int size = 0;
	stream >> size;
	assert(size > 1);
	try {
		matrix_t A(size, std::vector<double>(size, 0.0));
		read_matrix_from_file(stream, A);
		stream.close();
#ifdef DEBUG_PRINT
		PRINT_MATRIX(A);
#endif
		column_t f(size, 1);
		column_t x = find_solution(A, f);
		PRINT_COLUMN(x);

	}catch(const std::exception& err) {
		std::cerr << err.what();
		std::exit(EXIT_FAILURE);
	}
	return 0;
}

column_t find_solution(const matrix_t& A, const column_t& f)
{
	column_t x(f.size(), 0.0);
	column_t cur = x;
	double r = euclidean_norm_residual(A, f, x);
	std::cout << "Zero norm of residual " << r << std::endl;
	while (r >= ERROR) {

		for(int i = 0; i < f.size(); ++i){
			if (A[i][i] == 0)
				throw std::runtime_error("Division by zero!");
			cur[i] = (right_part(A, f, x, i) - previous_sum(A, cur, i)) / A[i][i];
		}

		x = cur;
		r = euclidean_norm_residual(A, f, x);

		std::cout << "Norm residual == " << r << std::endl;
	}

	return x;
}

double right_part(const matrix_t& A, const column_t& f, const column_t& x, int row)
{
	assert(row >= 0);
	double sum = f[row] * W - A[row][row] * (W - 1) * x[row];
	for (int i = row + 1; i < f.size(); ++i) {
		sum -= A[row][i] * W * x[i];
	}
	return sum;
}

double previous_sum(const matrix_t& A, const column_t& cur, int row)
{
	assert(row >= 0);
	double sum = 0.0;
	for (int i = 0; i < row; ++i) {
		sum += W * A[row][i] * cur[i];
	}
	return sum;
}

double euclidean_norm_residual(const matrix_t& A, const column_t& f, const column_t& x)
{
	column_t residual(f.size(), 0.0);
	double sum;
	for(int i = 0; i < f.size(); ++i) {
		sum = 0.0;
		for(int j = 0; j < f.size(); ++j) {
			sum += A[i][j] * x[j];
		}
		residual[i] = sum - f[i];
	}
	sum = 0.0;
	for(int i = 0; i < f.size(); ++i) {
		sum += (residual[i] * residual[i]);
	}
	return std::sqrt(sum);
}

void read_matrix_from_file(std::ifstream& in, matrix_t& M)
{
	int n = M.size();
	assert(n > 1);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			in >> M[i][j];
		}
	}
}