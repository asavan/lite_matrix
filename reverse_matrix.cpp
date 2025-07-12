#include "matrix.h"
#include "reverse.h"
#include <algorithm>
#include <math.h>

namespace matrix {
void Gauss(cmatr<double> A, cmatr<double>& B) {
    int size = A.get_size_stl();
    if (size != A.get_size_str()) {
        throw std::exception("can't use gauss method");
    }

    for (int i = 0; i < size; i++) {
        { // major element choise block
            double temp_max = 0;
            int addr = i;
            for (int k = i; k < size; k++) {
                if (fabs((double)A[k][i]) > fabs((double)temp_max)) {
                    temp_max = A[k][i];
                    addr = k;
                }
            }
            if (fabs((double)temp_max) < reverse::eps) {
                throw std::exception("det = 0");
            }
            if (addr != i) {
                // data_t *temp = A[i];	A[i] = A[addr];	A[addr] = temp;
                std::swap(A[i], A[addr]);
                // temp = B[i]; B[i] = B[addr]; B[addr] = temp;
                std::swap(B[i], B[addr]);
            }
        }
        for (int j = i + 1; j < size; j++) {
            double mn = A[j][i] / A[i][i];
            for (int k = i + 1; k < size; k++) {
                A[j][k] -= mn * A[i][k];
            }
            A[j][i] = 0;

            for (int k = 0; k < B.get_size_stl(); k++)
                B[j][k] -= mn * B[i][k];
        }
    }

    ///////////////////////////////////////reverse
    for (int i = size - 1; i >= 0; i--) {

        double mn = 1 / A[i][i];
        for (int k = 0; k < B.get_size_stl(); k++)
            B[i][k] *= mn;
        A[i][i] = 1;

        for (int j = i - 1; j >= 0; j--) {
            double mn1 = A[j][i];
            for (int k = 0; k < B.get_size_stl(); k++) {
                B[j][k] -= B[i][k] * mn1;
            }
            // for(int k = 0; k < size ; k++) {A[j][k]-=A[i][k]*mn1;} // chek for good work	}
        }
    }
}
cmatr<double> Gauss_Reverse(const cmatr<double>& A) {
    int size = A.get_size_stl();
    if (size != A.get_size_str()) {
        throw std::exception("can't use gauss method");
    }
    cmatr<double> B(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j)
                B[i][j] = 1;
            else
                B[i][j] = 0;
        }
    }
    Gauss(A, B);
    return B;
}
} // namespace matrix