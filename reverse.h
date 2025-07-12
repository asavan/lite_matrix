#ifndef __REVERSE_H__
#define __REVERSE_H__
#include "matrix.h"
namespace matrix {
cmatr<double> Gauss_Reverse(const cmatr<double>& A);
namespace reverse {
extern const double eps;
}
} // namespace matrix
#endif
