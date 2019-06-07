#include "timestepping.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
using namespace std::complex_literals;

template<typename T>
std::string str_real(T x) {
    std::ostringstream s;
    s << x;
    return s.str();
}
template<typename T>
std::string str_complex(std::complex<T> z) {
    std::ostringstream s;
    s << z.real() << " " << z.imag();
    return s.str();
}

template<typename T>
std::vector<T> deriv_exp(const std::vector<T>& y, T r=1) {
    std::vector<T> dy;
    dy.reserve(y.size());
    for(auto yi: y) {
        dy.push_back(r*yi);
    }
    return dy;
}

std::vector<double> deriv_oscillator(const std::vector<double>& y, double k=1) {
    return std::vector<double>{y[1], -k*y[0]};
}

template<typename dtype, typename stringer,
    typename DerivFn, typename Timestepper>
void write_ode(DerivFn deriv, const std::vector<dtype>& y0, double t_final,
    Timestepper step, std::string fname, stringer str) {
    
    std::ofstream outfile(fname);

    auto odeout = timestepping::odesolve(deriv, y0, t_final, step);
    for(auto point: odeout) {
        outfile << point.first;
        for(auto cmp: point.second) {
            outfile << " " << str(cmp);
        }
        outfile << std::endl;
    }
}

int main() {
    // auto step = timestepping::RK2(1e-3);
    // auto step = timestepping::RK4(1e-3);
    auto step = timestepping::AdaptiveRK(1e-3, 1e-1);

    double r = -1, k = 5, f = 1;
    double t_final = 10;
    write_ode([r](double t, const std::vector<double>& y) -> std::vector<double>
        {return deriv_exp(y, r);}, std::vector<double>{1, 0},
            t_final, step, "rexp.out", str_real<double>);
    
    write_ode([k](double t, const std::vector<double>& y) -> std::vector<double>
        {return deriv_oscillator(y, k);}, std::vector<double>{1, 0},
            t_final, step, "oscillator.out", str_real<double>);
    
    write_ode([f](double t, const std::vector<std::complex<double>>& y)
        -> std::vector<std::complex<double>> {
            return deriv_exp(y, 2.*M_PI*1i*f);
        },
        std::vector<std::complex<double>>{1},
            t_final, step, "cexp.out", str_complex<double>);
    write_ode([f](float t, const std::vector<std::complex<float>>& y)
        -> std::vector<std::complex<float>> {
            return deriv_exp(y, static_cast<std::complex<float>>(2.*M_PI*1i*f));
        },
        std::vector<std::complex<float>>{1},
            t_final, step, "cexpf.out", str_complex<float>);
}
