//
// Created by raj on 4/20/26.
//

#ifndef QISC_GATELIB_H
#define QISC_GATELIB_H

#include "core/gate.h"

/**
 * SINGLE QUBIT STATIC GATES
 */
Gate *CREATE_PAULI_X(void);      // Not-gate/Bit-flip
Gate *CREATE_PAULI_Y(void);
Gate *CREATE_PAULI_Z(void);      // Phase-flip
Gate *CREATE_HADAMARD(void);     // Superposition creator
Gate *CREATE_IDENTITY(void);     // The "Do Nothing" gate (actually useful for padding)
Gate *CREATE_S(void);            // Phase gate (sqrt of Z)
Gate *CREATE_S_DAGGER(void);     // Inverse Phase gate
Gate *CREATE_T(void);            // pi/8 gate
Gate *CREATE_T_DAGGER(void);     // Inverse pi/8 gate

/**
 * SINGLE QUBIT PARAMETRIC GATES
 * These require an angle (double theta)
 */
Gate *CREATE_RX(double theta);
Gate *CREATE_RY(double theta);
Gate *CREATE_RZ(double theta);
Gate *CREATE_PHASE(double phi);

/**
 * UNIVERSAL U-GATES (IBM Q Style)
 */
Gate *CREATE_U1(double lambda);
Gate *CREATE_U2(double phi, double lambda);
Gate *CREATE_U3(double theta, double phi, double lambda);

/**
 * TWO-QUBIT GATES
 */
Gate *CREATE_CNOT(void);         // Controlled-NOT (CX)
Gate *CREATE_CZ(void);           // Controlled-Z
Gate *CREATE_CPHASE(double theta);
Gate *CREATE_SWAP(void);         // Qubit swap
Gate *CREATE_ISWAP(void);        // Imaginary swap
Gate *CREATE_XX(void);           // Ising coupling gate
Gate *CREATE_YY(void);
Gate *CREATE_ZZ(void);

/**
 * THREE-QUBIT GATES
 */
Gate *CREATE_TOFFOLI(void);      // CCNOT (Universal for classical logic)
Gate *CREATE_FREDKIN(void);      // Controlled-SWAP





#endif //QISC_GATELIB_H
