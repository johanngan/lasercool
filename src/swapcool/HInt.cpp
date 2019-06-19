#include "HInt.hpp"
using namespace std::complex_literals;

unsigned HInt::subidx(unsigned i, unsigned j) const {
    return nstates*i + j;
}

std::vector<std::complex<double>> HInt::density_matrix(
    double gt, const std::vector<std::complex<double>>& coefficients) const {
    std::complex<double> cexp = std::exp(1i*cumulative_phase(gt));
    return {
        coefficients[subidx(0,0)],
        coefficients[subidx(0,1)],
        coefficients[subidx(0,2)],
        coefficients[subidx(1,0)],
        coefficients[subidx(1,1)],
        coefficients[subidx(1,2)]*cexp,
        coefficients[subidx(2,0)],
        coefficients[subidx(2,1)]*std::conj(cexp),
        coefficients[subidx(2,2)]
    };
}

// Assumes row major format
std::vector<std::complex<double>> HInt::operator()(double gt,
    const std::vector<std::complex<double>>& rho_c) const {
    // 1/(i*HBAR) * [H, rho_c] + L(rho_c) from the master equation
    return {
        (1 - branching_ratio)*rho_c[subidx(2,2)] * enable_decay,
        0.5i*(detun_per_decay(gt)*rho_c[subidx(0,1)]
            + rabi_softswitch(gt)*rho_c[subidx(0,2)]),
        -0.5*rho_c[subidx(0,2)] * enable_decay
            + 0.5i*(-detun_per_decay(gt)*rho_c[subidx(0,2)]
            + rabi_softswitch(gt)*rho_c[subidx(0,1)]),

        -0.5i*(detun_per_decay(gt)*rho_c[subidx(1,0)]
            + rabi_softswitch(gt)*rho_c[subidx(2,0)]),
        branching_ratio*rho_c[subidx(2,2)] * enable_decay
            + 0.5i*rabi_softswitch(gt)*(rho_c[subidx(1,2)]-rho_c[subidx(2,1)]),
        -0.5*rho_c[subidx(1,2)] * enable_decay
            + 1i*(0.5*rabi_softswitch(gt)*(rho_c[subidx(1,1)]-rho_c[subidx(2,2)])
            - detun_per_decay(gt)*rho_c[subidx(1,2)]),
        
        -0.5*rho_c[subidx(2,0)] * enable_decay
            + 0.5i*(detun_per_decay(gt)*rho_c[subidx(2,0)]
            - rabi_softswitch(gt)*rho_c[subidx(1,0)]),
        -0.5*rho_c[subidx(2,1)] * enable_decay
            - 1i*(0.5*rabi_softswitch(gt)*(rho_c[subidx(1,1)]-rho_c[subidx(2,2)])
            - detun_per_decay(gt)*rho_c[subidx(2,1)]),
        -rho_c[subidx(2,2)] * enable_decay
            - 0.5i*rabi_softswitch(gt)*(rho_c[subidx(1,2)]-rho_c[subidx(2,1)])
    };
}
