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
    double complex norm_sq = tensor_inner_product(circuit, circuit);
    if (!isfinite(creal(norm_sq)) || !isfinite(cimag(norm_sq))) {
        setlasterror("ERRGRV007: NON-FINITE STATE AMPLITUDE");
        return -1;
    }

    return 0;
}

static int validate_full_gate_capacity(const circuitt *circuit) {
    if (circuit->dim > SIZE_MAX / circuit->dim) {
        setlasterror("ERRGRV014: GROVER GATE MATRIX TOO LARGE");
        return -1;
    }

    size_t matrix_len = circuit->dim * circuit->dim;
    if (matrix_len > SIZE_MAX / sizeof(double complex)) {
        setlasterror("ERRGRV014: GROVER GATE MATRIX TOO LARGE");
        return -1;
    }

    return 0;
}

static int *create_all_targets(const circuitt *circuit) {
    int *targets = (int *)malloc((size_t)circuit->n * sizeof(int));
    if (!targets) {
        setlasterror("ERRGRV008: TARGET MEMORY ALLOCATION FAILED");
        return NULL;
    }

    for (int qubit = 0; qubit < circuit->n; qubit++) {
        targets[qubit] = qubit;
    }

    return targets;
}

static int apply_full_gate(circuitt *circuit, const Gate *gate) {
    int *targets = create_all_targets(circuit);
    if (!targets) {
        return -1;
    }

    int result = applygate(gate, circuit, targets);
    free(targets);
    if (result != 0) {
        return -1;
    }

    return validate_state_is_finite(circuit);
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

static Gate *create_zero_reflection_gate(const circuitt *circuit) {
    if (validate_full_gate_capacity(circuit) != 0) {
        return NULL;
    }

    Gate *gate = createGate(circuit->n);
    if (!gate) {
        return NULL;
    }

    for (size_t state = 1; state < gate->dim; state++) {
        gate->data[state * (gate->dim + 1U)] = -1.0;
    }

    return gate;
}

static Gate *create_phase_oracle_gate(const circuitt *circuit, grover_mark_fn oracle,
                                      void *user_data, size_t *out_marked_count) {
    if (validate_full_gate_capacity(circuit) != 0) {
        return NULL;
    }

    Gate *gate = createGate(circuit->n);
    if (!gate) {
        return NULL;
    }

    size_t marked_count = 0;
    for (size_t state = 0; state < circuit->dim; state++) {
        if (oracle(state, user_data)) {
            gate->data[state * (gate->dim + 1U)] = -1.0;
            marked_count++;
        }
    }

    if (marked_count == 0) {
        destroyGate(gate);
        setlasterror("ERRGRV006: ORACLE MARKED NO STATES");
        return NULL;
    }

    if (out_marked_count) {
        *out_marked_count = marked_count;
    }

    return gate;
}

static int apply_diffusion_reflection(circuitt *circuit, const Gate *reflection) {
    if (validate_state_is_finite(circuit) != 0) {
        return -1;
    }

    if (apply_all_hadamards(circuit) != 0) {
        return -1;
    }

    if (apply_full_gate(circuit, reflection) != 0) {
        return -1;
    }

    return apply_all_hadamards(circuit);
}

static int count_marked_states(const circuitt *circuit, grover_mark_fn oracle,
                               void *user_data, size_t *oco) {
    if (!oracle) {
        /// just in case
        setlasterror("ERRGRV005: NULL ORACLE");
        return -1;
    }

    size_t mcou = 0;
    for (size_t state = 0; state < circuit->dim; state++) {
        if (oracle(state, user_data)) {
            mcou++;
        }
    }

    if (mcou == 0) {
        setlasterror("ERRGRV006: ORACLE MARKED NO STATES");
        return -1;
    }

    if (oco) {
        *oco = mcou;
    }
    return 0;
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

    Gate *oracle_gate = create_phase_oracle_gate(circuit, oracle, user_data, NULL);
    if (!oracle_gate) {
        return -1;
    }

    int result = apply_full_gate(circuit, oracle_gate);
    destroyGate(oracle_gate);
    return result;
}

int groverdiffuse(circuitt *circuit) {
    if (validate_circuit(circuit) != 0) {
        return -1;
    }

    Gate *reflection = create_zero_reflection_gate(circuit);
    if (!reflection) {
        return -1;
    }

    int result = apply_diffusion_reflection(circuit, reflection);
    destroyGate(reflection);
    return result;
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
    if (iterations < 0) {
        size_t marked_count = 0;
        if (count_marked_states(circuit, oracle, user_data, &marked_count) != 0) {
            return -1;
        }

        resolved_iterations = groverrecit(circuit->n, marked_count);
        if (resolved_iterations < 0) {
            return -1;
        }
    }

    if (groveruniform(circuit) != 0) {
        return -1;
    }

    Gate *oracle_gate = NULL;
    Gate *reflection = NULL;
    if (resolved_iterations > 0) {
        oracle_gate = create_phase_oracle_gate(circuit, oracle, user_data, NULL);
        if (!oracle_gate) {
            return -1;
        }

        reflection = create_zero_reflection_gate(circuit);
        if (!reflection) {
            destroyGate(oracle_gate);
            return -1;
        }
    }

    for (int i = 0; i < resolved_iterations; i++) {
        if (apply_full_gate(circuit, oracle_gate) != 0) {
            destroyGate(oracle_gate);
            destroyGate(reflection);
            return -1;
        }
        if (apply_diffusion_reflection(circuit, reflection) != 0) {
            destroyGate(oracle_gate);
            destroyGate(reflection);
            return -1;
        }
    }

    destroyGate(oracle_gate);
    destroyGate(reflection);
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
