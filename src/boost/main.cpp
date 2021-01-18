#include <iostream>
#include <ctime>
#include <chrono>

#include <boost/numeric/odeint.hpp>

#include "swave.hpp"

using namespace boost::numeric::odeint;

struct push_back_delta_and_time {
    swave_quench& swq_;
    std::vector<double>& delta_;
    std::vector<double>& times_;

    push_back_delta_and_time(swave_quench &swq, std::vector<double> &delta, std::vector<double> &times) : swq_(swq), delta_(delta), times_(times) { };

    void operator() (const std::vector<double> &s, double t) {
        double d = swq_.gf_ * swq_.calc_delta(s);
        delta_.push_back(d);
        times_.push_back(t);
    }
};

int main() {

    double D0 = 0.6;
    double Df = 0.8;
    double D = 10.0;
    int Nspins = 50000;

    swave_quench swq(D0, Df, D, Nspins);

    double t0 = 0.0;
    double tf = 100.0;
    int num_time_steps = 10000;
    double dt = tf / num_time_steps;
    double reltol = 1.0e-4;
    double abstol = 1.0e-9;
    int max_num_steps = 10000;
    size_t NEQ = 3 * Nspins;

    std::vector<double> s0(NEQ);
    swq.set_initial_state(s0);

    std::vector<double> delta;
    std::vector<double> times;

    printf("#N: %d\n", swq.Nspins_);
    printf("#D0: %.12f\n", swq.D0_);
    printf("#Df: %.12f\n", swq.Df_);
    printf("#g0: %.12f\n", swq.g0_);
    printf("#gf: %.12f\n", swq.gf_);

    typedef runge_kutta_dopri5<std::vector<double>> error_stepper_type;

    auto start_time = std::chrono::steady_clock::now();
    //size_t steps = integrate(swq, s0, t0, tf, dt, push_back_delta_and_time(swq, delta, times));
    //size_t steps = integrate_times(make_dense_output<error_stepper_type>(abstol, reltol), swq, s0, t0, tf, dt, push_back_delta_and_time(swq, delta,times));
    //size_t steps = integrate_adaptive(make_controlled<error_stepper_type>(abstol, reltol), swq, s0, t0, tf, dt, push_back_delta_and_time(swq, delta,times));
    size_t steps = integrate_adaptive(make_dense_output<error_stepper_type>(abstol, reltol), swq, s0, t0, tf, dt, push_back_delta_and_time(swq, delta,times));
    auto end_time = std::chrono::steady_clock::now();

    std::cout << "#Elapsed time: " << std::chrono::duration<double>(end_time - start_time).count() << " seconds\n";

    for (size_t i = 0; i <= steps; i++) {
        std::cout << times[i] << "," << delta[i] << "\n";
    }

    return 0;
}
