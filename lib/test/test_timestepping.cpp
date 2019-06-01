#include "timestepping.hpp"
#include <fstream>
#include <string>
#include <vector>

std::vector<double> deriv_growth(std::vector<double>, double r=1);
std::vector<double> deriv_oscillator(std::vector<double>, double k=1);

template<typename DerivFn, typename Timestepper>
void write_ode(DerivFn deriv, std::vector<double> y0, double t_final,
    Timestepper step, std::string fname) {
    
    std::ofstream outfile(fname);

    auto odeout = timestepping::odesolve(deriv, y0, t_final, step);
    for(auto point: odeout) {
        outfile << point.first;
        for(auto cmp: point.second) {
            outfile << " " << cmp;
        }
        outfile << std::endl;
    }
}

int main() {
    auto step = timestepping::AdaptiveRK(1e-6, 1e-1);

    double r = -1, k = 5;
    double t_final = 10;
    write_ode([r](double t, std::vector<double> y) -> std::vector<double>
        {return deriv_growth(y, r);}, {1, 0}, t_final, step, "growth.out");
    
    write_ode([k](double t, std::vector<double> y) -> std::vector<double>
        {return deriv_oscillator(y, k);}, {1, 0}, t_final, step, "oscillator.out");
}

std::vector<double> deriv_growth(std::vector<double> y, double r) {
    std::vector<double> dy;
    dy.reserve(y.size());
    for(auto yi: y) {
        dy.push_back(r*yi);
    }
    return dy;
}

std::vector<double> deriv_oscillator(std::vector<double> y, double k) {
    return std::vector<double>{y[1], -k*y[0]};
}