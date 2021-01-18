#ifndef SWAVE_HPP
#define SWAVE_HPP

#include <vector>

class swave_quench {

    public:
    double D0_;       // Initial coupling
    double g0_;       // Initial coupling
    double Df_;       // Final coupling
    double gf_;       // Final coupling
    double D_;        // Half bandwidth
    int Nspins_;      // Number of spins
    std::vector<double> epsk_;    // energy levels

    swave_quench(double D0, double Df, double D, int Nspins) : D0_(D0), Df_(Df), D_(D), Nspins_(Nspins) {
        for (int i = 0; i < Nspins; i++) {
            epsk_.push_back(-D + 2 * D * i / (Nspins - 1));
        }
        g0_ = calc_g(D0);
        gf_ = calc_g(Df);
    }
    void set_initial_state(std::vector<double> &s);

    double calc_g(const double delta);
    double calc_delta(const std::vector<double> &s);
    //void eom(const std::vector<double> &y, std::vector<double> &ydot, const double t)
    void operator() (const std::vector<double> &y, std::vector<double> &ydot, const double t);

};

#endif
