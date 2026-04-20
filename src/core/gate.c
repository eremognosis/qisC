//
// Created by raj on 4/20/26.
//
#include "core/gate.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "core/error.h"

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
    gate->data = data;
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
