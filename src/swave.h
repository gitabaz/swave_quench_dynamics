#ifndef SWAVE_H
#define SWAVE_H

#include <cvode/cvode.h>

typedef struct {
    double D0;       // Initial coupling
    double g0;       // Initial coupling
    double Df;       // Final coupling
    double gf;       // Final coupling
    double D;        // Half bandwidth
    int Nspins;      // Number of spins
    double* epsk;    // energy levels
} quench_params;

quench_params* swave_quench_init(double D0, double Df, double D, int Nspins);
void swave_set_initial_state(double *s, quench_params *qp);
void swave_set_epsk(quench_params *qp);
void swave_free_quench_params(quench_params *qp);

double swave_calc_g(char c, quench_params* qp);
double swave_calc_delta(double *s, quench_params *qp);
int swave_eom(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#endif
