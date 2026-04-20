#ifndef QISC_GROVER_H
#define QISC_GROVER_H

#include "core/qubit.h"
#include <stdbool.h>
#include <stddef.h>

typedef bool (*grover_mark_fn)(size_t state, void *user_data);

int groverrecit(int n_qubits, size_t marked_count);
int groveruniform(circuitt *circuit);
int grover_apply_oracle(circuitt *circuit, grover_mark_fn oracle, void *user_data);
int groverdiffuse(circuitt *circuit);
int grover_run(circuitt *circuit, grover_mark_fn oracle, void *user_data, int iterations);
int grover_search(circuitt *circuit, size_t marked_state, int iterations);
int grovermlstate(const circuitt *circuit, size_t *out_state, double *out_probability);

#endif
