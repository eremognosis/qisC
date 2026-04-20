//
// Created by raj on 4/20/26.
//

#ifndef QISC_GATE_H
#define QISC_GATE_H
#include <complex.h>
#include <stddef.h>
#include "core/qubit.h"
typedef struct {
    int n;
    size_t dim; /// 2**n but just redundancy soi that a osychic aint eqquired to underfstand
    double complex *data; /// flattenned
} Gate;



Gate *createGate(int n); //// oresents to I
void destroyGate(Gate *gate);
int setGate(Gate *gate, double complex *data);
int applygate(const Gate *gate, circuitt *circuit, const int *target_qubits); /// wtf


#endif //QISC_GATE_H
