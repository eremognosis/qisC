#ifndef QISC_TENSOR_H
#define QISC_TENSOR_H

#include "core/qubit.h"
#include <complex.h>
#include <stddef.h>




int tensor_kronecker(const double complex *A, size_t rA, size_t cA,
                            const double complex *B, size_t rB, size_t cB,
                            double complex **out_C);

double complex tensor_inner_product(const circuitt *psi, const circuitt *phi);
double complex tensor_trace(const double complex *mat, size_t dim);


int tensor_partial_trace(const double complex *rho_in, int n_qubits,
                                int target_qubit, double complex *rho_out);

void tensor_normalize(circuitt *circuit);

#endif //QISC_TENSOR_H