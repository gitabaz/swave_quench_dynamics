#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>

//#include <cvode/cvode.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_butcher_erk.h>
#include <nvector/nvector_openmp.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "swave.h"

int main() {

    double D0 = 0.6;
    double Df = 0.8;
    double D = 10.0;
    int Nspins = 50000;
    quench_params *qp = swave_quench_init(D0, Df, D, Nspins);

    double t0 = 0.0;
    double tf = 100.0;
    int num_time_steps = 10000;
    double reltol = 1.0e-4;
    double abstol = 1.0e-9;
    int max_num_steps = 10000;

    double *delta = malloc((num_time_steps + 1) * sizeof(*delta));

    sunindextype NEQ = 3 * Nspins;

    N_Vector s0 = N_VNew_OpenMP(NEQ, 4);
    realtype *s0data = NV_DATA_OMP(s0);

    swave_set_initial_state(s0data, qp);

    //void *cvode_mem = CVodeCreate(CV_ADAMS);
    void *arkode_mem = ERKStepCreate(swave_eom, t0, s0);
    
    int flag = 0;
    //flag = CVodeInit(cvode_mem, swave_eom, t0, s0);
    //flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    flag = ERKStepSStolerances(arkode_mem, reltol, abstol);

    //CVodeSetUserData(cvode_mem, (void *) qp);
    ERKStepSetUserData(arkode_mem, (void *) qp);
    ERKStepSetTableNum(arkode_mem, DORMAND_PRINCE_7_4_5);

    //SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(s0, 0);
    //flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
    //flag = CVodeSetMaxNumSteps(cvode_mem, max_num_steps);

    double tout = 0.0;
    double dt = tf / num_time_steps;

    printf("#N: %d\n", qp->Nspins);
    printf("#D0: %.12f\n", qp->D0);
    printf("#Df: %.12f\n", qp->Df);
    printf("#g0: %.12f\n", qp->g0);
    printf("#gf: %.12f\n", qp->gf);

    time_t start_time = time(NULL);
    delta[0] = qp->gf * swave_calc_delta(s0data, qp);
    for (int i = 1; i <= num_time_steps; i++) {
        //flag = CVode(cvode_mem, i * dt, s0, &tout, CV_NORMAL);
        flag = ERKStepEvolve(arkode_mem, i * dt, s0, &tout, ARK_NORMAL);
        delta[i] = qp->gf * swave_calc_delta(s0data, qp);
    }
    time_t end_time = time(NULL);

    printf("#Elapsed time: %jd seconds\n", end_time - start_time);

    for (int i = 0; i <= num_time_steps; i++) {
        printf("%.4f,%.12f\n", t0 + i * dt, delta[i]);
    }


    free(delta);
    swave_free_quench_params(qp);
    N_VDestroy(s0);
    ERKStepFree(&arkode_mem);
    //CVodeFree(&cvode_mem);
    //SUNNonlinSolFree(NLS);
    
    return 0;
}
