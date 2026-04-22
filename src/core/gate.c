//
// Created by raj on 4/20/26.
//
#include "core/gate.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "core/error.h"

#define GATE_UNITARY_EPSILON 1e-10
#define GATE_ORTHO_EPSILON 1e-14

static int gate_matrix_len(const Gate *gate, size_t *out_len) {
    if (!gate) {
        setlasterror("ERRGT004: NULL GATE");
        return -1;
    }
    if (gate->n < 1 || gate->n >= (int)(sizeof(size_t) * CHAR_BIT)) {
        setlasterror("ERRGT009: INVALID GATE QUBIT COUNT");
        return -1;
    }

    size_t expected_dim = ((size_t)1U) << gate->n;
    if (gate->dim != expected_dim) {
        setlasterror("ERRGT012: INVALID GATE DIMENSION");
        return -1;
    }
    if (gate->dim > SIZE_MAX / gate->dim ||
        gate->dim * gate->dim > SIZE_MAX / sizeof(double complex)) {
        setlasterror("ERRGT017: GATE MATRIX TOO LARGE");
        return -1;
    }

    if (out_len) {
        *out_len = gate->dim * gate->dim;
    }
    return 0;
}

static double complex column_dot(const double complex *matrix, size_t dim,
                                 size_t left_col, const double complex *right_vec) {
    double complex dot = 0.0;
    for (size_t row = 0; row < dim; row++) {
        dot += conj(matrix[row * dim + left_col]) * right_vec[row];
    }
    return dot;
}

static double vector_norm(const double complex *vec, size_t dim) {
    double norm_sq = 0.0;
    for (size_t row = 0; row < dim; row++) {
        double real = creal(vec[row]);
        double imag = cimag(vec[row]);
        if (!isfinite(real) || !isfinite(imag)) {
            return NAN;
        }
        norm_sq += real * real + imag * imag;
    }
    return sqrt(norm_sq);
}

static bool matrix_is_finite(const double complex *matrix, size_t len) {
    for (size_t i = 0; i < len; i++) {
        if (!isfinite(creal(matrix[i])) || !isfinite(cimag(matrix[i]))) {
            return false;
        }
    }
    return true;
}

static void remove_previous_columns(const double complex *unitary, size_t dim,
                                    size_t upto_col, double complex *vec) {
    for (size_t prev_col = 0; prev_col < upto_col; prev_col++) {
        double complex projection = column_dot(unitary, dim, prev_col, vec);
        for (size_t row = 0; row < dim; row++) {
            vec[row] -= unitary[row * dim + prev_col] * projection;
        }
    }
}

static void transpose(double complex *mat, size_t dim) {
    /// a matrix of size dim * dim flattened
    if (dim<2) {
        return;
    }
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = i + 1; j < dim; j++) {
            size_t idx1 = i * dim + j;
            size_t idx2 = j * dim + i;

            double complex temp = mat[idx1];
            mat[idx1] = mat[idx2];
            mat[idx2] = temp;
        }
    }
}

Gate *createGate(int n) {
    if (n == 0) {
        setlasterror("ERRGT001: GATE WITH 0 DIMENSION");
        return NULL;
    }
    Gate *gate = (Gate *) calloc(1, sizeof(Gate));
    if (!gate) {
        setlasterror("ERRGT002: MEMORY ALLOCATION FAILED");
        return NULL;
    }
    gate->n = n;
    gate->dim = ((size_t)1U)<<n;
    gate->data = (double complex *)calloc(gate->dim * gate->dim, sizeof(double complex));

    if (!gate->data) {
        setlasterror("ERRGT003: DATA MEMORY ALLOCATION FAILED");
        free(gate);
        return NULL;
    }
    for (size_t i = 0; i < gate->dim; i++) {
        // 0, dim+1, 2*dim+2,
        gate->data[i*(gate->dim+1)] = 1.00;

    }
    return gate;
}

void destroyGate(Gate *gate) {
    free(gate->data);
    free(gate);
}

int setGate(Gate *gate, double complex *data) {
    if (!gate) {
        setlasterror("ERRGT004: NULL GATE");
        return -1;
    }
    if (!data) {
        setlasterror("ERRGT005: NULL DATA");
        return -1;
    }

    size_t matrix_len = 0;
    if (gate_matrix_len(gate, &matrix_len) != 0) return -1;

    if (!matrix_is_finite(data, matrix_len)) {
        setlasterror("ERRGT019: NON-FINITE GATE MATRIX VALUE");
        return -1;
    }

    /// we create a ghost matrix actually to see if shit works
    Gate ghost_gate = *gate;
    ghost_gate.data = data;

    if (!checkunitary(&ghost_gate, GATE_UNITARY_EPSILON)) {
        setlasterror("ERRGT021: GATE IS NOT UNITARY");
        return -1; /////// original gate remains untouched and safe.
    }

    ///// that we know the probability (and other things) are not leaking we actually commit to memory (well we have to commit some things to memory)
    if (!gate->data) {
        gate->data = (double complex *)malloc(matrix_len * sizeof(double complex));
        if (!gate->data) {
            setlasterror("ERRGT003: DATA MEMORY ALLOCATION FAILED");
            return -1;
        }
    }

    memcpy(gate->data, data, matrix_len * sizeof(double complex));
    return 0;
}

int applygate(const Gate *gate, circuitt *circuit, const int *target_qubits) {
    if (!gate) {
        setlasterror("ERRGT004: NULL GATE");
        return -1;
    }
    if (!circuit) {
        setlasterror("ERRGT005: NULL CIRCUIT");

        return -1;
    }
    if (!target_qubits) {
        setlasterror("ERRGT006: NULL TARGET QUBITS");
        return -1;
    }
    if (!gate->data) {
        setlasterror("ERRGT007: NULL GATE DATA");
        return -1;
    }
    if (!circuit->state_amps) {
        setlasterror("ERRGT008: NULL CIRCUIT STATE");
        return -1;
    }
    if (gate->n < 1 || gate->n >= (int)(sizeof(size_t) * CHAR_BIT)) {
        setlasterror("ERRGT009: INVALID GATE QUBIT COUNT");
        return -1;
    }
    if (circuit->n < 1 || circuit->n >= (int)(sizeof(size_t) * CHAR_BIT)) {
        setlasterror("ERRGT010: INVALID CIRCUIT QUBIT COUNT");
        return -1;
    }
    if (gate->n > circuit->n) {
        /// this aint gfonna happen like the mmeory will fuck much faster
        setlasterror("ERRGT011: GATE LARGER THAN CIRCUIT");
        return -1;
    }

    size_t egdim = ((size_t)1U) << gate->n;
    size_t ecdim = ((size_t)1U) << circuit->n;
    if (gate->dim != egdim) {
        setlasterror("ERRGT012: INVALID GATE DIMENSION");
        return -1;
    }
    if (circuit->dim != ecdim) {
        setlasterror("ERRGT013: INVALID CIRCUIT DIMENSION");
        return -1;
    }

    size_t target_mask = 0;
    for (int i = 0; i < gate->n; i++) {
        if (target_qubits[i] < 0 || target_qubits[i] >= circuit->n) {
            setlasterror("ERRGT014: TARGET QUBIT OUT OF RANGE");
            return -1;
        }

        size_t bit = ((size_t)1U) << target_qubits[i];
        if (target_mask & bit) {
            setlasterror("ERRGT015: DUPLICATE TARGET QUBIT");
            return -1;
        }
        target_mask |= bit;
    }

    if (gate->n == 1) {
        const double complex m00 = gate->data[0];
        const double complex m01 = gate->data[1];
        const double complex m10 = gate->data[2];
        const double complex m11 = gate->data[3];

        for (size_t base = 0; base < circuit->dim;) {
            size_t i0 = base;
            size_t i1 = base | target_mask;
            double complex a0 = circuit->state_amps[i0];
            double complex a1 = circuit->state_amps[i1];

            circuit->state_amps[i0] = m00 * a0 + m01 * a1;
            circuit->state_amps[i1] = m10 * a0 + m11 * a1;

            size_t next = (base | target_mask) + 1U;
            if (next >= circuit->dim) {
                break;
            }
            base = next & ~target_mask;
        }

        return 0;
    }

    if (gate->dim > SIZE_MAX / sizeof(size_t) ||
        gate->dim > SIZE_MAX / sizeof(double complex)) {
        setlasterror("ERRGT016: STATE MEMORY ALLOCATION FAILED");
        return -1;
    }

    size_t *basis_offsets = (size_t *)malloc(gate->dim * sizeof(size_t));
    double complex *input = (double complex *)malloc(gate->dim * sizeof(double complex));
    double complex *output = (double complex *)malloc(gate->dim * sizeof(double complex));
    if (!basis_offsets || !input || !output) {
        free(basis_offsets);
        free(input);
        free(output);
        setlasterror("ERRGT016: STATE MEMORY ALLOCATION FAILED");
        return -1;
    }

    for (size_t gate_basis = 0; gate_basis < gate->dim; gate_basis++) {
        size_t circuit_offset = 0;
        for (int target = 0; target < gate->n; target++) {
            if (gate_basis & (((size_t)1U) << target)) {
                circuit_offset |= ((size_t)1U) << target_qubits[target];
            }
        }
        basis_offsets[gate_basis] = circuit_offset;
    }

    for (size_t base = 0; base < circuit->dim;) {
        for (size_t col = 0; col < gate->dim; col++) {
            input[col] = circuit->state_amps[base | basis_offsets[col]];
        }

        for (size_t row = 0; row < gate->dim; row++) {
            double complex amp = 0.0;
            for (size_t col = 0; col < gate->dim; col++) {
                amp += gate->data[row * gate->dim + col] * input[col];
            }
            output[row] = amp;
        }

        for (size_t row = 0; row < gate->dim; row++) {
            circuit->state_amps[base | basis_offsets[row]] = output[row];
        }

        size_t next = (base | target_mask) + 1U;
        if (next >= circuit->dim) {
            break;
        }
        base = next & ~target_mask;
    }

    free(basis_offsets);
    free(input);
    free(output);
    return 0;
}


bool checkunitary(const Gate *gate, double epsilon) {
    if (gate_matrix_len(gate, NULL) != 0) {
        return false;
    }
    if (!gate->data) {
        setlasterror("ERRGT007: NULL GATE DATA");
        return false;
    }
    if (!isfinite(epsilon) || epsilon < 0.0) {
        setlasterror("ERRGT018: INVALID UNITARY CHECK EPSILON");
        return false;
    }

    for (size_t left_col = 0; left_col < gate->dim; left_col++) {
        for (size_t right_col = 0; right_col < gate->dim; right_col++) {
            double complex dot = 0.0;
            for (size_t row = 0; row < gate->dim; row++) {
                dot += conj(gate->data[row * gate->dim + left_col]) *
                       gate->data[row * gate->dim + right_col];
            }

            if (!isfinite(creal(dot)) || !isfinite(cimag(dot))) {
                setlasterror("ERRGT019: NON-FINITE GATE MATRIX VALUE");
                return false;
            }

            double complex expected = (left_col == right_col) ? 1.0 : 0.0;
            if (cabs(dot - expected) > epsilon) {
                return false;
            }
        }
    }

    return true;
}

int makeunitary(Gate *gate) {
    size_t matrix_len = 0;
    if (gate_matrix_len(gate, &matrix_len) != 0) {
        return -1;
    }
    if (!gate->data) {
        setlasterror("ERRGT007: NULL GATE DATA");
        return -1;
    }
    if (!matrix_is_finite(gate->data, matrix_len)) {
        setlasterror("ERRGT019: NON-FINITE GATE MATRIX VALUE");
        return -1;
    }
    if (checkunitary(gate, GATE_UNITARY_EPSILON)) {
        return 0;
    }

    double complex *unitary = (double complex *)calloc(matrix_len, sizeof(double complex));
    double complex *vec = (double complex *)malloc(gate->dim * sizeof(double complex));
    if (!unitary || !vec) {
        free(unitary);
        free(vec);
        setlasterror("ERRGT016: STATE MEMORY ALLOCATION FAILED");
        return -1;
    }

    for (size_t col = 0; col < gate->dim; col++) {
        for (size_t row = 0; row < gate->dim; row++) {
            vec[row] = gate->data[row * gate->dim + col];
        }

        remove_previous_columns(unitary, gate->dim, col, vec);
        double norm = vector_norm(vec, gate->dim);

        if (!isfinite(norm) || norm <= GATE_ORTHO_EPSILON) {
            bool found_basis = false;

            for (size_t basis_row = 0; basis_row < gate->dim; basis_row++) {
                for (size_t row = 0; row < gate->dim; row++) {
                    vec[row] = (row == basis_row) ? 1.0 : 0.0;
                }

                remove_previous_columns(unitary, gate->dim, col, vec);
                norm = vector_norm(vec, gate->dim);
                if (isfinite(norm) && norm > GATE_ORTHO_EPSILON) {
                    found_basis = true;
                    break;
                }
            }

            if (!found_basis) {
                free(unitary);
                free(vec);
                setlasterror("ERRGT020: COULD NOT MAKE GATE UNITARY");
                return -1;
            }
        }

        for (size_t row = 0; row < gate->dim; row++) {
            unitary[row * gate->dim + col] = vec[row] / norm;
        }
    }

    memcpy(gate->data, unitary, matrix_len * sizeof(double complex));
    free(unitary);
    free(vec);

    if (!checkunitary(gate, GATE_UNITARY_EPSILON)) {
        setlasterror("ERRGT020: COULD NOT MAKE GATE UNITARY");
        return -1;
    }

    return 0;
}

Gate *gatedagger(Gate *gate) {
    if (!gate) {
        setlasterror("ERRGT004: NULL GATE DATA");
        return NULL;
    }
    if (!gate->data) {
        setlasterror("ERRGT005: NULL GATE DATA");
        return NULL;
    }
    int n = gate->n;
    size_t dim = gate->dim;
    Gate *newgate= createGate(n);
    if (!newgate) {
        setlasterror("ERRGT0023: GATE ALLOCATION FAILED");
        return NULL;
    }

    setGate(newgate, gate->data);

    for (size_t i = 0; i < dim*dim; i++) {
        newgate->data[i] = conj(newgate->data[i]);
    }

    transpose(newgate->data, newgate->dim);

    return newgate;
}

