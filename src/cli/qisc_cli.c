#include "qisc.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(FILE *stream, const char *program) {
    fprintf(stream,
            "Usage:\n"
            "  %s bell\n"
            "  %s grover <n_qubits> <marked_state> [iterations]\n"
            "  %s --help\n"
            "\n"
            "Examples:\n"
            "  %s bell\n"
            "  %s grover 3 5\n",
            program, program, program, program, program);
}

static int parse_int_arg(const char *text, const char *name, int *out_value) {
    char *end = NULL;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        fprintf(stderr, "Invalid %s: %s\n", name, text);
        return -1;
    }

    *out_value = (int)value;
    return 0;
}

static int parse_size_arg(const char *text, const char *name, size_t *out_value) {
    char *end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        fprintf(stderr, "Invalid %s: %s\n", name, text);
        return -1;
    }

    *out_value = (size_t)value;
    return 0;
}

static void print_last_error(void) {
    const char *message = getlasterror();
    if (message && message[0] != '\0') {
        fprintf(stderr, "%s\n", message);
    }
}

static int print_statevector(const circuitt *circuit) {
    if (!circuit) {
        return -1;
    }

    qScT *state = (qScT *)calloc(circuit->dim, sizeof(qScT));
    if (!state) {
        fprintf(stderr, "Could not allocate state-vector output buffer.\n");
        return -1;
    }

    if (statevec(circuit, state, circuit->dim) != 0) {
        print_last_error();
        free(state);
        return -1;
    }

    for (size_t basis = 0; basis < circuit->dim; basis++) {
        double probability = state[basis].real * state[basis].real +
                             state[basis].imag * state[basis].imag;
        printf("|%zu>  amplitude=% .12f%+.12fi  probability=%.12f\n",
               basis, state[basis].real, state[basis].imag, probability);
    }

    free(state);
    return 0;
}

static int run_bell(void) {
    circuitt *circuit = createcirc(2);
    Gate *hadamard = CREATE_HADAMARD();
    Gate *cnot = CREATE_CNOT();
    if (!circuit || !hadamard || !cnot) {
        print_last_error();
        destroycirc(circuit);
        destroyGate(hadamard);
        destroyGate(cnot);
        return 1;
    }

    int h_target = 0;
    int cnot_targets[2] = {0, 1};
    if (applygate(hadamard, circuit, &h_target) != 0 ||
        applygate(cnot, circuit, cnot_targets) != 0) {
        print_last_error();
        destroycirc(circuit);
        destroyGate(hadamard);
        destroyGate(cnot);
        return 1;
    }

    printf("Bell circuit: H(0), CNOT(0, 1)\n");
    int result = print_statevector(circuit) == 0 ? 0 : 1;

    destroycirc(circuit);
    destroyGate(hadamard);
    destroyGate(cnot);
    return result;
}

static int run_grover(int argc, char **argv) {
    if (argc < 4 || argc > 5) {
        print_usage(stderr, argv[0]);
        return 1;
    }

    int n_qubits = 0;
    size_t marked_state = 0;
    int iterations = -1;
    if (parse_int_arg(argv[2], "n_qubits", &n_qubits) != 0 ||
        parse_size_arg(argv[3], "marked_state", &marked_state) != 0) {
        return 1;
    }
    if (argc == 5 && parse_int_arg(argv[4], "iterations", &iterations) != 0) {
        return 1;
    }

    int resolved_iterations = iterations;
    if (resolved_iterations < 0) {
        resolved_iterations = groverrecit(n_qubits, 1);
        if (resolved_iterations < 0) {
            print_last_error();
            return 1;
        }
    }

    circuitt *circuit = createcirc(n_qubits);
    if (!circuit) {
        print_last_error();
        return 1;
    }

    if (grover_search(circuit, marked_state, iterations) != 0) {
        print_last_error();
        destroycirc(circuit);
        return 1;
    }

    size_t measured_state = 0;
    double probability = 0.0;
    if (grovermlstate(circuit, &measured_state, &probability) != 0) {
        print_last_error();
        destroycirc(circuit);
        return 1;
    }

    printf("Grover search\n");
    printf("  qubits: %d\n", n_qubits);
    printf("  marked state: %zu\n", marked_state);
    printf("  iterations: %d\n", resolved_iterations);
    printf("  most likely state: %zu\n", measured_state);
    printf("  probability: %.12f\n", probability);

    if (circuit->dim <= 256) {
        printf("\nState vector:\n");
        if (print_statevector(circuit) != 0) {
            destroycirc(circuit);
            return 1;
        }
    } else {
        printf("\nState vector omitted because it has %zu amplitudes.\n", circuit->dim);
    }

    destroycirc(circuit);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        print_usage(argc < 2 ? stderr : stdout, argv[0]);
        return argc < 2 ? 1 : 0;
    }

    if (strcmp(argv[1], "bell") == 0) {
        return run_bell();
    }
    if (strcmp(argv[1], "grover") == 0) {
        return run_grover(argc, argv);
    }

    fprintf(stderr, "Unknown command: %s\n", argv[1]);
    print_usage(stderr, argv[0]);
    return 1;
}
