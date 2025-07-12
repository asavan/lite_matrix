#include "matrix.h"
#include "reverse.h"
#include <iostream>

const double eps = 1e-14;
int matrix::default_matrix_size = 2;
typedef matrix::cmatr<double> matr;
	
int main(int argc, char * argv[])
{
	using namespace matrix;	
	matr A("matrix.txt");
	std::cout<<Gauss_Reverse(A);
	return 0;
}
