#include <cmath>

#include "swave.hpp"

void swave_quench::set_initial_state(std::vector<double> &s) {
    double Ek;
    for (int i = 0; i < Nspins_; i++) {
        Ek = std::sqrt(std::pow(epsk_[i], 2) + std::pow(D0_, 2));
        s[i] = D0_ / (2.0 * Ek);                                 // Sx
        s[i + Nspins_] = 0.0;                                    // Sy
        s[i + 2 * Nspins_] = -epsk_[i] / (2.0 * Ek);             // Sz
    }
}

double swave_quench::calc_delta(const std::vector<double> &s) {
    double delta = 0.0;
    for (int i = 0; i < Nspins_; i++) {
        delta += s[i];
    }
    return delta;
}

double swave_quench::calc_g(const double delta) {
    double res = 0.0;
    for (int i = 0; i < Nspins_; i++) {
        res += 1.0 / std::sqrt(std::pow(epsk_[i], 2) + std::pow(delta, 2));
    }
    return 2.0 / res;
}

void swave_quench::operator() (const std::vector<double> &y, std::vector<double> &ydot, const double t) {
    double delta = gf_ * calc_delta(y);
#pragma omp parallel for
    for (int i = 0; i < Nspins_; i++) {
        ydot[i]                  = - 2 * epsk_[i] * y[i + Nspins_];
        ydot[i + Nspins_]     =   2 * epsk_[i] * y[i] + 2 * delta * y[i + 2 * Nspins_];
        ydot[i + 2 * Nspins_] = - 2 * delta * y[i + Nspins_];
    }
}
