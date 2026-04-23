//
// Created by raj on 4/21/26.
//
#include "algos/grover.h"

#include "core/error.h"
#include "core/tensor.h"
#include "gates/gatelib.h"
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#define GROVER_PI 3.141592653589793238462643383279502884

typedef struct {
    size_t marked_state;
} single_mark_ctx;

static bool mark_single_state(size_t state, void *user_data) {
    single_mark_ctx *ctx = (single_mark_ctx *)user_data;
    return ctx && state == ctx->marked_state;
}

static int validate_circuit(const circuitt *circuit) {
    if (!circuit) {
        setlasterror("ERRGRV001: NULL CIRCUIT");
        return -1;
    }
    if (!circuit->state_amps) {
        setlasterror("ERRGRV002: NULL CIRCUIT STATE");
        return -1;
    }
    if (circuit->n < 1 || circuit->n >= (int)(sizeof(size_t) * CHAR_BIT)) {
        setlasterror("ERRGRV003: INVALID CIRCUIT QUBIT COUNT");
        return -1;
    }

    size_t expected_dim = ((size_t)1U) << circuit->n;
    if (circuit->dim != expected_dim) {
        setlasterror("ERRGRV004: INVALID CIRCUIT DIMENSION");
        return -1;
    }

    return 0;
}

static int validate_state_is_finite(const circuitt *circuit) {
    for (size_t state = 0; state < circuit->dim; state++) {
        double real = creal(circuit->state_amps[state]);
        double imag = cimag(circuit->state_amps[state]);
        if (!isfinite(real) || !isfinite(imag)) {
            setlasterror("ERRGRV007: NON-FINITE STATE AMPLITUDE");
            return -1;
        }
    }

    return 0;
}

static int gather_marked_states(const circuitt *circuit, grover_mark_fn oracle,
                                void *user_data, size_t **out_states,
                                size_t *out_count) {
    if (!oracle) {
        setlasterror("ERRGRV005: NULL ORACLE");
        return -1;
    }
    if (!out_states || !out_count) {
        setlasterror("ERRGRV013: NULL OUTPUT STATE");
        return -1;
    }

    size_t *marked_states = NULL;
    size_t marked_count = 0;
    size_t marked_capacity = 0;

    for (size_t state = 0; state < circuit->dim; state++) {
        if (!oracle(state, user_data)) {
            continue;
        }

        if (marked_count == marked_capacity) {
            size_t new_capacity = (marked_capacity == 0) ? 8U : marked_capacity * 2U;
            if (new_capacity < marked_capacity ||
                new_capacity > SIZE_MAX / sizeof(size_t)) {
                free(marked_states);
                setlasterror("ERRGRV008: TARGET MEMORY ALLOCATION FAILED");
                return -1;
            }

            size_t *resized = (size_t *)realloc(marked_states,
                                                new_capacity * sizeof(size_t));
            if (!resized) {
                free(marked_states);
                setlasterror("ERRGRV008: TARGET MEMORY ALLOCATION FAILED");
                return -1;
            }

            marked_states = resized;
            marked_capacity = new_capacity;
        }

        marked_states[marked_count++] = state;
    }

    if (marked_count == 0) {
        free(marked_states);
        setlasterror("ERRGRV006: ORACLE MARKED NO STATES");
        return -1;
    }

    *out_states = marked_states;
    *out_count = marked_count;
    return 0;
}

static void apply_marked_phase_flip(circuitt *circuit, const size_t *marked_states,
                                    size_t marked_count) {
    for (size_t i = 0; i < marked_count; i++) {
        size_t marked_state = marked_states[i];
        circuit->state_amps[marked_state] = -circuit->state_amps[marked_state];
    }
}

static int apply_diffusion_reflection(circuitt *circuit) {
    if (validate_state_is_finite(circuit) != 0) {
        return -1;
    }

    double complex average = 0.0;
    for (size_t state = 0; state < circuit->dim; state++) {
        average += circuit->state_amps[state];
    }
    average /= (double)circuit->dim;

    for (size_t state = 0; state < circuit->dim; state++) {
        circuit->state_amps[state] = (2.0 * average) - circuit->state_amps[state];
    }

    return validate_state_is_finite(circuit);
}

static int apply_oracle_phase_flip(circuitt *circuit, grover_mark_fn oracle,
                                   void *user_data) {
    size_t *marked_states = NULL;
    size_t marked_count = 0;
    if (gather_marked_states(circuit, oracle, user_data, &marked_states,
                             &marked_count) != 0) {
        return -1;
    }

    apply_marked_phase_flip(circuit, marked_states, marked_count);
    free(marked_states);

    if (validate_state_is_finite(circuit) != 0) {
        return -1;
    }

    return 0;
}

static int apply_all_hadamards(circuitt *circuit) {
    Gate *hadamard = CREATE_HADAMARD();
    if (!hadamard) {
        return -1;
    }

    for (int qubit = 0; qubit < circuit->n; qubit++) {
        int target = qubit;
        if (applygate(hadamard, circuit, &target) != 0) {
            destroyGate(hadamard);
            return -1;
        }
    }

    destroyGate(hadamard);
    return validate_state_is_finite(circuit);
}

int groverrecit(int n_qubits, size_t marked_count) {
    if (n_qubits < 1 || n_qubits >= (int)(sizeof(size_t) * CHAR_BIT)) {
        setlasterror("ERRGRV009: INVALID QUBIT COUNT");
        return -1;
    }

    size_t dim = ((size_t)1U) << n_qubits;
    if (marked_count == 0 || marked_count > dim) {
        setlasterror("ERRGRV010: INVALID MARKED STATE COUNT");
        return -1;
    }
    if (marked_count == dim) {
        return 0;
    }

    double marked_ratio = (double)marked_count / (double)dim;
    double theta = asin(sqrt(marked_ratio));
    double estimate = (GROVER_PI / (4.0 * theta)) - 0.5;
    if (estimate <= 0.0) {
        return 0;
    }
    if (estimate > (double)INT_MAX) {
        setlasterror("ERRGRV011: RECOMMENDED ITERATION COUNT TOO LARGE");
        return -1;
    }

    return (int)lround(estimate);
}

int groveruniform(circuitt *circuit) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }

    if (resetcirc(circuit) != 0) {
        return -1;
    }
    circuit->state_amps[0] = 1.0;

    if (apply_all_hadamards(circuit) != 0) {
        return -1;
    }

    tensor_normalize(circuit);
    return validate_state_is_finite(circuit);
}

int grover_apply_oracle(circuitt *circuit, grover_mark_fn oracle, void *user_data) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }
    if (!oracle) {
        setlasterror("ERRGRV005: NULL ORACLE");
        return -1;
    }

    return apply_oracle_phase_flip(circuit, oracle, user_data);
}

int groverdiffuse(circuitt *circuit) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }

    return apply_diffusion_reflection(circuit);
}

int grover_run(circuitt *circuit, grover_mark_fn oracle, void *user_data, int iterations) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }
    if (!oracle) {
        setlasterror("ERRGRV005: NULL ORACLE");
        return -1;
    }

    int resolved_iterations = iterations;
    size_t *marked_states = NULL;
    size_t marked_count = 0;
    if (iterations < 0) {
        if (gather_marked_states(circuit, oracle, user_data, &marked_states,
                                 &marked_count) != 0) {
            return -1;
        }

        resolved_iterations = groverrecit(circuit->n, marked_count);
        if (resolved_iterations < 0) {
            free(marked_states);
            return -1;
        }
    }
    else if (iterations > 0 &&
             gather_marked_states(circuit, oracle, user_data, &marked_states,
                                  &marked_count) != 0) {
        return -1;
    }

    if (groveruniform(circuit) != 0) {
        free(marked_states);
        return -1;
    }

    for (int i = 0; i < resolved_iterations; i++) {
        apply_marked_phase_flip(circuit, marked_states, marked_count);
        if (validate_state_is_finite(circuit) != 0) {
            free(marked_states);
            return -1;
        }
        if (apply_diffusion_reflection(circuit) != 0) {
            free(marked_states);
            return -1;
        }
    }

    free(marked_states);
    tensor_normalize(circuit);
    return validate_state_is_finite(circuit);
}

int grover_search(circuitt *circuit, size_t marked_state, int iterations) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }
    if (marked_state >= circuit->dim) {
        setlasterror("ERRGRV012: MARKED STATE OUT OF RANGE");
        return -1;
    }

    single_mark_ctx ctx = {marked_state};
    return grover_run(circuit, mark_single_state, &ctx, iterations);
}

int grovermlstate(const circuitt *circuit, size_t *out_state, double *out_probability) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }
    if (!out_state) {
        setlasterror("ERRGRV013: NULL OUTPUT STATE");
        return -1;
    }
    if (validate_state_is_finite(circuit) != 0) {
        return -1;
    }

    size_t best_state = 0;
    double best_probability = -1.0;
    for (size_t state = 0; state < circuit->dim; state++) {
        double magnitude = cabs(circuit->state_amps[state]);
        double probability = magnitude * magnitude;
        if (probability > best_probability) {
            best_state = state;
            best_probability = probability;
        }
    }

    *out_state = best_state;
    if (out_probability) {
        *out_probability = best_probability;
    }

    return 0;
}
