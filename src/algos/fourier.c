//
// Created by raj on 4/22/26.
//

#include "algos/fourier.h"
#include "core/gate.h"
#include "core/error.h"
#include "core/qubit.h"
#include "gates/gatelib.h"
#include <math.h>
#include <stdlib.h>
#define FOURIER_PI 3.141592653589793238462643383279502884
// The quantum Fourier transform is the classical discrete Fourier transform applied to the vector of amplitudes of a quantum state, which has length {\displaystyle N=2^{n}} if it is applied to a register of  {\displaystyle n} qubits



void applyQFT(circuitt* circuit, int *target, int nq) {
    if (!circuit) { setlasterror("ERRFRR001: NULL CIRCUIT"); return; }
    if (!circuit->state_amps) { setlasterror("ERRFRR002: EMOTY CIRCUIT"); return; }
    if (!target) { setlasterror("ERRFRR003: NULL target"); return; }
    if (nq > circuit->n) { setlasterror("ERRFR004: target more than SIZE"); return; }

    const Gate *hada = CREATE_HADAMARD();

    for (int i = 0; i < nq; i++) {
        //  Hadamard to the current qubit
        applygate(hada, circuit, &target[i]);

        for (int j = i + 1; j < nq; j++) {
            double theta = (2.0 * FOURIER_PI) / pow(2, j - i + 1);
            Gate *cph = CREATE_CPHASE(theta);
            int qubits[2] = {target[j], target[i]};
            applygate(cph, circuit, qubits);
            destroyGate(cph);
        }
    }

    const Gate *swp = CREATE_SWAP();
    for (int i = 0; i < nq / 2; i++) {
        int swap_targets[2] = {target[i], target[nq - i - 1]};
        applygate(swp, circuit, swap_targets);
    }

    destroyGate(swp);
    destroyGate(hada);
}


void applyiQFT(circuitt* circuit, int *target, int nq) {

    const Gate *swp = CREATE_SWAP();
    for (int i = nq / 2 - 1; i >= 0; i--) {
        int swap_targets[2] = {target[i], target[nq - i - 1]};
        applygate(swp, circuit, swap_targets);
    }
    destroyGate(swp);
    const Gate *hada = CREATE_HADAMARD();
    for (int i = nq - 1; i >= 0; i--) {

        for (int j = nq - 1; j > i; j--) {

            double theta = -(2.0 * FOURIER_PI) / (double)(1 << (j - i + 1));
            Gate *ch = CREATE_CPHASE(theta);
            int qubits[2] = {target[j], target[i]};
            applygate(ch, circuit, qubits);
            destroyGate(ch);
        }


        applygate(hada,circuit,&target[i]);

    }
    destroyGate(hada);
}