#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>
/*
 * The task is to solve non-linear system of equations of such a type
 * u = A * e^(-u) via Newton Method
 * at each iteration we are to solve a system of linear algebraic equations via Gauss method
 * */
using column_t = std::vector<double>;
using matrix_t = std::vector<column_t>;

/*
 * the base condition for an iterative method to stop
 * calculating new approximation of the desirable solution
 * */
static const double ERROR = 1e-6;

matrix_t read_matrix_from_file(std::ifstream& in);
void print_matrix(const matrix_t &A);
column_t gauss_method(matrix_t& A, column_t& b);
double sum_previous_elems(const matrix_t &A, const column_t &x, unsigned int i);
double euclidean_norm(const column_t &x);
void residual(const matrix_t &A, const column_t &x0, column_t &r);
void compute_jacobi(matrix_t &jacobi, const column_t &x0, const matrix_t &A);
void solution(const matrix_t &A);

int main(int argc, char **argv)
{
	if (argc != 2) {
		std::cout << "Bad amount of input arguments\n";
		std::exit(EXIT_FAILURE);
	}

	std::ifstream in(argv[1]);
	if (!in.is_open()) {
		std::cout << "Error in opening file " << argv[1] << std::endl;
		std::exit(EXIT_FAILURE);
	}
	matrix_t A = read_matrix_from_file(in);
	in.close();
	try{
		solution(A);
	}catch(const std::exception& err) {
		std::cout << err.what();
		std::exit(EXIT_FAILURE);
	}
	return 0;
}

void solution(const matrix_t &A)
{
	column_t x0(A.size(), 0); // zero approximation
	//f(u) = u - A * e^(-u)
	matrix_t jacobi = A; // we just initialize it with some values
	column_t delta = x0;
	column_t r = x0;
	residual(A, x0, r);
	double norm = euclidean_norm(r);
	std::cout << "Zero norm residual == " << norm << std::endl;

	int iter = 0;
	while (norm >= ERROR) {
		compute_jacobi(jacobi, x0, A);
		for (auto &el : r)
			el = -el;
		delta = gauss_method(jacobi, r); // delta = x_(k+1) - x_k
		for (unsigned int i = 0; i < delta.size(); ++i)
			x0[i] += delta[i];
		residual(A, x0, r);
		norm = euclidean_norm(r);
		std::cout << "[ " << iter << " ] " << "Norm residual: " << norm << std::endl;
		std::cout << "[ " << iter << " ] " << "Current approximation: ";
		for (const auto &item : x0)
			std::cout << item << " ";
		std::cout << std::endl;
		++iter;
	}
}

void residual(const matrix_t &A, const column_t &x0, column_t &r)
{
	//compute a vector: f(x0) = x0 - A * e^(-x0)
	int id = 0, i = 0;
	r = x0;
	for (const auto &line : A) {
		for (const auto &el : line) {
			r[id] -= el * std::exp(-x0[i]);
			++i;
		}
		i = 0;
		++id;
	}
}

void compute_jacobi(matrix_t &jacobi, const column_t &x0, const matrix_t &A)
{
	for(unsigned int i = 0; i < jacobi.size(); ++i) {
		for (unsigned int j = 0; j < jacobi.size(); ++j) {
			if (i == j)
				jacobi[i][i] = 1 + A[i][i] * std::exp(-x0[i]);
			else
				jacobi[i][j] = A[i][j] * std::exp(-x0[j]);
		}
	}
}

double euclidean_norm(const column_t &x)
{
	double res = 0;
	for (const auto &el : x) {
		res += el * el;
	}
	return std::sqrt(res);
}

matrix_t read_matrix_from_file(std::ifstream& in)
{
	assert(in.good());
	matrix_t A;
	int size = 0;
	in >> size;
	assert(size > 1);
	A = matrix_t(size, column_t(size, 0.0)); // move assignment operator
	for (auto &line : A) {
		for (auto &el : line)
			in >> el;
	}
	return A;
}

void print_matrix(const matrix_t& A)
{
	for (const auto &line : A) {
		for (const auto &el : line)
			std::cout << el << " ";
		std::cout << std::endl;
	}
}

column_t gauss_method(matrix_t& A, column_t& b)
{
	double coeff = 0;
	for (unsigned int line = 0; line < (A.size() - 1); ++line) {
		if (A[line][line] == 0) { // here we suppose that the first element in a row isn't null
			throw std::runtime_error("Division by zero");
		}
		for (unsigned int line_below = (line + 1); line_below < A.size(); ++line_below) {
			if (A[line_below][line] == 0)
				continue;
			coeff = A[line_below][line] / A[line][line];
			for (unsigned int el = line; el < A.size(); ++el) {
				A[line_below][el] -= coeff * A[line][el];

			}
			b[line_below] -= b[line] * coeff;
			//assert(A[line_below][line] == 0);
		}
	}
	//reverse gauss method
	column_t x(A.size(), 0.0);
	for (int i = static_cast<int>(A.size() - 1); i >= 0; --i) {
		if (A[i][i] == 0)
			throw std::runtime_error("Division by zero");
		x[i] = (b[i] - sum_previous_elems(A, x, i)) / A[i][i];
	}
	return x;
}

double sum_previous_elems(const matrix_t &A, const column_t &x, unsigned int i)
{
	double res = 0;
	for (unsigned int j = i + 1; j < A.size(); ++j) {
		res += x[j] * A[i][j];
	}
	return res;
}