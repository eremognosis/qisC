#include "core/tensor.h"
#include "core/error.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>



double complex tensor_inner_product(const circuitt *psi, const circuitt *phi) {
    // vcalcucaltes the <psi | phi >

    if (!psi || !phi) {
        setlasterror("ERRTNSR001: STATEVECTOR IS NULL");
        return CMPLX(NAN, NAN);
    }

    if (!psi->state_amps || !phi->state_amps) {
        setlasterror("ERRTNSR002: STATEVECTOR IS EMPTY");
        return CMPLX(NAN, NAN);
    }


    size_t dim1 = psi->dim;
    size_t dim = phi->dim;
    if (dim1 != dim) {
        setlasterror("ERRTNSR003: STATE VECTOR OF DIFFERENT DIMENSIONS");
        return CMPLX(NAN, NAN);
    }
    double complex res = 0.0;
    for (size_t i = 0; i < dim; i++) { //why the fuck index sdhit looks okay at 4 am and breaks at 4:58
        res += conj(psi->state_amps[i] )* phi->state_amps[i];
    }

    return res;

}


double complex tensor_trace(const double complex *mat, size_t dim) {
    if (!mat || !dim) {
        setlasterror("ERRTNSR004: INVALID INPUT");
        return CMPLX(NAN, NAN);
    }

    double complex res = 0.0;
    for (size_t i = 0; i < dim; i++) {
        res += *(mat+i+i*dim);
    }

    return res;

}


int tensor_kronecker(const double complex *A, size_t rA, size_t cA,
                     const double complex *B, size_t rB, size_t cB,
                     double complex **out_C) {
    if (!A || !B || !out_C || !rA || !cA || !rB || !cB) {
        setlasterror("ERRTNSR004: INVALID INPUT");
        return -1;
    }

    for (size_t r = 0; r < rA; r++) {
        for (size_t c = 0; c < cA; c++) {
            double complex av = A[r * cA + c];

            for (size_t ib = 0; ib < rB; ib++) {
                for (size_t jb = 0; jb < cB; jb++) {
                    size_t rC = r * rB + ib;
                    size_t cC = c * cB + jb;

                    out_C[rC][cC] = av * B[ib * cB + jb];
                }
            }
        }
    }

    return 0;
}


int tensor_partial_trace(const double complex *rho_in, int n_qubits,
                         int target_qubit, double complex *rho_out) {
    if (!rho_out || !rho_in) {
        setlasterror("ERRTNSR004: INVALID INPUT");
        return -1;
    }


    size_t news = ((size_t)1U) << (n_qubits - 1);

    size_t olds = ((size_t)1U) << n_qubits;

    size_t low_mask = ((size_t)1U << target_qubit) - 1;

    for (size_t i = 0; i < news; i++) {
        size_t r0 = ((i & ~low_mask) << 1) | (0 << target_qubit) | (i & low_mask);
        size_t r1 = r0 | (1 << target_qubit);

        for (size_t j = 0; j < news; j++) {
            size_t c0 = ((j & ~low_mask) << 1) | (0 << target_qubit) | (j & low_mask);
            size_t c1 = c0 | (1 << target_qubit);
            rho_out[i * news + j] = rho_in[r0 * olds + c0] + rho_in[r1 * olds + c1];
        }
    }

    return 0;
}

void tensor_normalize(circuitt *circuit) {
    if (!circuit) {
        setlasterror("ERRTNSR001: STATEVECTOR IS NULL");
        return;
    }

    if (!circuit->state_amps) {
        setlasterror("ERRTNSR002: STATEVECTOR IS EMPTY");
        return;
    }

    double complex A = tensor_inner_product(circuit, circuit);  //normalization coinstant

    if (!A) {
        setlasterror("ERRTNSR006: ZERO VECTOR");
        return;
    }

    size_t dim = circuit->dim;
    double complex norm = csqrt(A);
    for (size_t i = 0; i < dim; i++) {
        circuit->state_amps[i] /= norm;
    }

}

