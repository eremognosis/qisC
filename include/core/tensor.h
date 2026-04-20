#ifndef QISC_TENSOR_H
#define QISC_TENSOR_H

#include "core/qubit.h"
#include <complex.h>
#include <stddef.h>




typedef struct {
    double complex data[2][2];
} Gate2x2;


int tensor_kronecker(const double complex *A, size_t rA, size_t cA,
                            const double complex *B, size_t rB, size_t cB,
                            double complex **out_C);

double complex tensor_inner_product(const circuitt *psi, const circuitt *phi);
double complex tensor_trace(const double complex *mat, size_t dim);

/**
 * @brief Partial Trace: rho_A = Tr_B(rho_AB)
 * Reduces a density matrix by tracing out the specified qubit.
 */
int tensor_partial_trace(const double complex *rho_in, int n_qubits,
                                int target_qubit, double complex *rho_out);

/**
 * @brief Normalizes the state vector to ensure sum(|amp|^2) == 1.
 * Essential to call periodically to fight floating point drift.
 */
void tensor_normalize(circuitt *circuit);

#endif //QISC_TENSOR_H