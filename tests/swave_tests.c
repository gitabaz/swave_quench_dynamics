#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

#include "../src/swave.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#define GREEN "\033[92m"
#define NOCOLOR "\033[0m"

static int tests_run = 0;

static void test_quench_init() {
    tests_run++;
    double D0 = 0.5;
    double Df = 1.0;
    double D = 10.0;
    int Nspins = 1000;

    quench_params *qp = swave_quench_init(D0, Df, D, Nspins);

    assert(qp->D0 == D0);
    assert(qp->Df == Df);
    assert(qp->D == D);
    assert(qp->Nspins == Nspins);

    swave_free_quench_params(qp);
}

static void test_swave_initial_state() {
    tests_run++;
    double D0 = 0.5;
    double Df = 1.0;
    double D = 10.0;
    int Nspins = 10000;
    quench_params *qp = swave_quench_init(D0, Df, D, Nspins);

    int neq = 3 * Nspins;
    double *s = malloc(neq * sizeof(*s));
    swave_set_initial_state(s, qp);

    double g0 = swave_calc_g(D0, qp);

    double Delta = g0 * swave_calc_delta(s, qp);

    assert(fabs(Delta - D0) < 1.0e-4);

    swave_free_quench_params(qp);
    free(s);
}

static int cvode_eom(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    /* van der pol oscillator
     * dy1 = y2
     * dy2 = 1000 * (1 - y1^2) * y2 - y1
     */
    realtype *ydata = N_VGetArrayPointer(y);
    realtype *ydotdata = N_VGetArrayPointer(ydot);

    ydotdata[0] = ydata[1];
    ydotdata[1] = 1000.0 * (1.0 - pow(ydata[0], 2)) * ydata[1] - ydata[0];

    return 0;
}

static void test_cvode() {
    tests_run++;
    sunindextype NEQ = 2;
    int max_num_steps = 10000;
    N_Vector y0 = N_VNew_Serial(NEQ);
    realtype *ydata = N_VGetArrayPointer(y0);
    ydata[0] = 2.0;
    ydata[1] = 0.0;

    //void* cvode_mem = CVodeCreate(CV_ADAMS);
    void* cvode_mem = CVodeCreate(CV_BDF);
    assert(cvode_mem);

    int flag = CVodeInit(cvode_mem, cvode_eom, 0.0, y0);
    assert(flag == CV_SUCCESS);

    realtype reltol = 1.0e-10;
    realtype abstol = 1.0e-10;

    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    assert(flag == CV_SUCCESS);

    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ);
    SUNLinearSolver LS = SUNLinSol_Dense(y0, A);
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    assert(flag == CV_SUCCESS);

    CVodeSetMaxNumSteps(cvode_mem, max_num_steps);
    
    realtype tf = 3000.0;
    realtype tret = 0.0;
    int num_steps = 10000;
    realtype dt = tf / num_steps;
    for (int i = 1; i <= num_steps; i++) {
        flag = CVode(cvode_mem, i * dt, y0, &tret, CV_NORMAL);
    }

    assert((ydata[0] - (-1.5106068)) < 1.0e-4);
    
    N_VDestroy(y0);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
}

static void test_calc_g() {
    tests_run++;
    double D0 = 0.5;
    double Df = 0.5;
    double D = 10.0;
    int Nspins = 50000;
    quench_params *qp = swave_quench_init(D0, Df, D, Nspins);

    double g0 = swave_calc_g(D0, qp);
    double gf = swave_calc_g(Df, qp);
    
    assert(fabs(g0 - 0.000114044380) < 1.0e-4);
    assert(fabs(gf - 0.000114044380) < 1.0e-4);

    swave_free_quench_params(qp);
}

static void test_swave_dynamics() {
    tests_run++;

    double D0 = 0.6;
    double Df = 0.8;
    double D = 10.0;
    int Nspins = 1000;
    quench_params *qp = swave_quench_init(D0, Df, D, Nspins);

    double t0 = 0.0;
    double tf = 10.0;
    int num_time_steps = 100;
    double reltol = 1.0e-10;
    double abstol = 1.0e-10;
    int max_num_steps = 10000;

    sunindextype NEQ = 3 * Nspins;

    N_Vector s0 = N_VNew_Serial(NEQ);
    realtype *s0data = N_VGetArrayPointer(s0);

    swave_set_initial_state(s0data, qp);

    void *cvode_mem = CVodeCreate(CV_BDF);

    int flag = 0;
    flag = CVodeInit(cvode_mem, swave_eom, t0, s0);
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);

    CVodeSetUserData(cvode_mem, (void *) qp);

    SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(s0, 1);
    flag = CVodeSetNonlinearSolver(cvode_mem, NLS);

    double tout = 0.0;
    double dt = tf / num_time_steps;
    double delta = 0.0;

    for (int i = 1; i <= num_time_steps; i++) {
        flag = CVode(cvode_mem, i * dt, s0, &tout, CV_NORMAL);
    }
    delta = qp->gf * swave_calc_delta(s0data, qp);

    assert(fabs(delta - 0.809964079314) < 1.0e-4);

    swave_free_quench_params(qp);
    N_VDestroy(s0);
    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
}

int main() {

    test_cvode();
    test_quench_init();
    test_swave_initial_state();
    test_calc_g();
    test_swave_dynamics();
    printf(GREEN "Passed %d tests\n" NOCOLOR, tests_run);

    return 0;
}
