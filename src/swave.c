#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

#include <nvector/nvector_openmp.h>

#include "swave.h"

quench_params* swave_quench_init(double D0, double Df, double D, int Nspins) {
    quench_params *qp = malloc(sizeof(*qp));

    if (qp) {
        qp->D0 = D0;
        qp->Df = Df;
        qp->D = D;
        qp->Nspins = Nspins;
        qp->epsk = malloc(Nspins * sizeof(*(qp->epsk)));
        swave_set_epsk(qp);
        qp->g0 = swave_calc_g(D0, qp);
        qp->gf = swave_calc_g(Df, qp);
    }

    return qp;
}

void swave_free_quench_params(quench_params *qp) {
    free(qp->epsk);
    free(qp);
}

void swave_set_epsk(quench_params *qp) {
    for (int i = 0; i < qp->Nspins; i++) {
        qp->epsk[i] = -qp->D + 2 * qp->D * i / (qp->Nspins - 1);
    }
}

void swave_set_initial_state(double *s, quench_params *qp) {
    double Ek;
    for (int i = 0; i < qp->Nspins; i++) {
        Ek = sqrt(pow(qp->epsk[i], 2) + pow(qp->D0, 2));
        s[i] = qp->D0 / (2.0 * Ek);                                 // Sx
        s[i + qp->Nspins] = 0.0;                                    // Sy
        s[i + 2 * qp->Nspins] = -qp->epsk[i] / (2.0 * Ek);          // Sz
    }
}

inline double swave_calc_delta(double* s, quench_params *qp) {
    double delta = 0.0;
    for (int i = 0; i < qp->Nspins; i++) {
        delta += s[i];
    }
    return delta;
}

double swave_calc_g(double delta, quench_params* qp) {
    double res = 0.0;
    for (int i = 0; i < qp->Nspins; i++) {
        res += 1.0 / sqrt(pow(qp->epsk[i], 2) + pow(delta, 2));
    }
    return 2.0 / res;
}

inline int swave_eom(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype *ydata = NV_DATA_OMP(y);
    realtype *ydotdata = NV_DATA_OMP(ydot);

    quench_params *qp = (quench_params *) user_data;

    double delta = qp->gf * swave_calc_delta(ydata, qp);
    for (int i = 0; i < qp->Nspins; i++) {
        ydotdata[i]                  = - 2 * qp->epsk[i] * ydata[i + qp->Nspins];
        ydotdata[i + qp->Nspins]     =   2 * qp->epsk[i] * ydata[i] + 2 * delta * ydata[i + 2 * qp->Nspins];
        ydotdata[i + 2 * qp->Nspins] = - 2 * delta * ydata[i + qp->Nspins];
    }
    return 0;
}
