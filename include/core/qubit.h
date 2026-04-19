#ifndef QSIM_QUBIT_H
#define QSIM_QUBIT_H

#include <complex.h>
#include <stddef.h>

typedef struct {
    double real;
    double imag;
} qScT;  /// complex

typedef struct {
    int n;
    size_t dim;
    double complex *state_amps;
} circuitt;

circuitt *createcirc(int n);
void destroycirc(circuitt *circuit);
int resetcirc(circuitt *circuit);
int statevec(const circuitt *circuit, qScT *out_state, size_t out_len);

#endif
