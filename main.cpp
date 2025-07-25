#include "matrix.h"
#include "reverse.h"
#include <iostream>

const double matrix::reverse::eps = 1e-14;
int matrix::default_matrix_size = 2;
typedef matrix::cmatr<double> matr;

int main(int /*argc*/, char* /*argv*/[]) {
    using namespace matrix;
    try {
        matr A("matrix.txt");
        std::cout << Gauss_Reverse(A);
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
