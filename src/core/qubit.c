//
// Created by raj on 4/20/26.
//
#include "core/qubit.h"
#include "core/error.h"
#include <string.h>
#include <stdlib.h>
#include <complex.h>
circuitt *createcirc(int n) {
    clearerror();

    if (n<1 || n>30) {
        setlasterror("ERRQBT001: n must be in 1 to 30, got %d\n", n); ///idk hoiw i will implentb error, lkets see later
        return NULL;
    }

    circuitt *circuit = (circuitt*) calloc(1, sizeof(circuitt));

    if (!circuit) {
        setlasterror("ERRQBT002: Memory Allocation Failed\n");
        return NULL;
    }

    circuit->n = n;
    circuit->dim = ((size_t)1U)<<n; //// very fancy wauy of saying 2**n
    circuit->state_amps = (double complex *)calloc(circuit->dim, sizeof(double complex));

    if (!circuit->state_amps) {
        free(circuit);
        setlasterror("ERRQBT003: Memory Allocation Failed for the state amplitudes\n");
        return NULL;
    }

    circuit->state_amps[0] = 1.00 + 0.00*I;

    return circuit;


}



void destroycirc(circuitt *circuit) {
    if (!circuit) {return;} ///silent reyturn works asd many times wre do redundant thingd
    free(circuit->state_amps);
    free(circuit);
}


int resetcirc(circuitt *circuit) {
    if (!circuit || !circuit->state_amps) {
        setlasterror("ERRQBT004: RESET CIRCUIT FAILED. NULL CIRCUIT");
        return -1;
    }

    memset(circuit->state_amps, 0, circuit->dim*sizeof(double complex));
    return 0;
}


int statevec(const circuitt *circuit, qScT *out_state, size_t out_len) {
    if (!circuit || !circuit->state_amps || !out_state) {
        setlasterror("ERRQBT005:");
        return -1;
    }
    if (out_len<circuit->dim) {
        setlasterror("ERRQBT006: STATEVECTOR EXPOPT BUFFER TOO SMALL. NEEDED %zu, got %zu",circuit->dim, out_len);
        return -1;
    }

    for (size_t i=0; i<circuit->dim; i++) {
        out_state[i].real = creal(circuit->state_amps[i]);
        out_state[i].imag = cimag(circuit->state_amps[i]);
    }
    return 0;
}

