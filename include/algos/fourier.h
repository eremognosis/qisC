//
// Created by raj on 4/22/26.
//

#ifndef QISC_FOURIER_H
#define QISC_FOURIER_H

#include "core/qubit.h"
#include "core/gate.h"


void applyQFT(circuitt* circuit, int *target, int nq);
void applyiQFT(circuitt* circuit, int *target, int nq);

#endif //QISC_FOURIER_H
