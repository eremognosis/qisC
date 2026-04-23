//
// Created by raj on 4/20/26.
//
#include "gates/gatelib.h"

#include <complex.h>
#include <math.h>

#define QISC_PI 3.141592653589793238462643383279502884

static double complex phase(double angle) {
    return cos(angle) + I * sin(angle);
}

static Gate *create_from_data(const int n, const double complex *data) {
    Gate *gate = createGate(n);
    if (!gate) {
        return NULL;
    }

    if (setGateTrusted(gate, data) != 0) {
        destroyGate(gate);
        return NULL;
    }

    return gate;
}

Gate *CREATE_PAULI_X(void) {
    double complex data[4] = {
        0.0, 1.0,
        1.0, 0.0,
    };
    return create_from_data(1, data);
}

Gate *CREATE_PAULI_Y(void) {
    double complex data[4] = {
        0.0, -I,
        I, 0.0,
    };
    return create_from_data(1, data);
}

Gate *CREATE_PAULI_Z(void) {
    double complex data[4] = {
        1.0, 0.0,
        0.0, -1.0,
    };
    return create_from_data(1, data);
}

Gate *CREATE_HADAMARD(void) {
    double s = 1.0 / sqrt(2.0);
    double complex data[4] = {
        s, s,
        s, -s,
    };
    return create_from_data(1, data);
}

Gate *CREATE_IDENTITY(void) {
    return createGate(1);
}

Gate *CREATE_S(void) {
    double complex data[4] = {
        1.0, 0.0,
        0.0, I,
    };
    return create_from_data(1, data);
}

Gate *CREATE_S_DAGGER(void) {
    double complex data[4] = {
        1.0, 0.0,
        0.0, -I,
    };
    return create_from_data(1, data);
}

Gate *CREATE_T(void) {
    double complex data[4] = {
        1.0, 0.0,
        0.0, phase(QISC_PI / 4.0),
    };
    return create_from_data(1, data);
}

Gate *CREATE_T_DAGGER(void) {
    double complex data[4] = {
        1.0, 0.0,
        0.0, phase(-QISC_PI / 4.0),
    };
    return create_from_data(1, data);
}

Gate *CREATE_RX(double theta) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    double complex data[4] = {
        c, -I * s,
        -I * s, c,
    };
    return create_from_data(1, data);
}

Gate *CREATE_RY(double theta) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    double complex data[4] = {
        c, -s,
        s, c,
    };
    return create_from_data(1, data);
}

Gate *CREATE_RZ(double theta) {
    double complex data[4] = {
        phase(-theta / 2.0), 0.0,
        0.0, phase(theta / 2.0),
    };
    return create_from_data(1, data);
}

Gate *CREATE_PHASE(double phi) {
    double complex data[4] = {
        1.0, 0.0,
        0.0, phase(phi),
    };
    return create_from_data(1, data);
}

Gate *CREATE_U1(double lambda) {
    return CREATE_PHASE(lambda);
}

Gate *CREATE_U2(double phi, double lambda) {
    double s = 1.0 / sqrt(2.0);
    double complex e_phi = phase(phi);
    double complex e_lambda = phase(lambda);
    double complex data[4] = {
        s, -e_lambda * s,
        e_phi * s, e_phi * e_lambda * s,
    };
    return create_from_data(1, data);
}

Gate *CREATE_U3(double theta, double phi, double lambda) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    double complex e_phi = phase(phi);
    double complex e_lambda = phase(lambda);
    double complex data[4] = {
        c, -e_lambda * s,
        e_phi * s, e_phi * e_lambda * c,
    };
    return create_from_data(1, data);
}

Gate *CREATE_CNOT(void) {
    double complex data[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_CZ(void) {
    double complex data[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, -1.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_CPHASE(double theta) {
    double complex data[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, phase(theta),
    };
    return create_from_data(2, data);
}

Gate *CREATE_SWAP(void) {
    double complex data[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_ISWAP(void) {
    double complex data[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, I, 0.0,
        0.0, I, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_XX(void) {
    double complex data[16] = {
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_YY(void) {
    double complex data[16] = {
        0.0, 0.0, 0.0, -1.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0, 0.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_ZZ(void) {
    double complex data[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, -1.0, 0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    };
    return create_from_data(2, data);
}

Gate *CREATE_TOFFOLI(void) {
    double complex data[64] = {0};
    for (size_t i = 0; i < 8; i++) {
        data[i * 8 + i] = 1.0;
    }
    data[3 * 8 + 3] = 0.0;
    data[7 * 8 + 7] = 0.0;
    data[7 * 8 + 3] = 1.0;
    data[3 * 8 + 7] = 1.0;
    return create_from_data(3, data);
}

Gate *CREATE_FREDKIN(void) {
    double complex data[64] = {0};
    for (size_t i = 0; i < 8; i++) {
        data[i * 8 + i] = 1.0;
    }
    data[3 * 8 + 3] = 0.0;
    data[5 * 8 + 5] = 0.0;
    data[5 * 8 + 3] = 1.0;
    data[3 * 8 + 5] = 1.0;
    return create_from_data(3, data);
}
