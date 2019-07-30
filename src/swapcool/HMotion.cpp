#include "HMotion.hpp"
using namespace std::complex_literals;

HMotion::HMotion(std::string fname):HSwap(fname), stationary_decay_prob(0.6) {}

std::vector<std::complex<double>> HMotion::density_matrix(
    double gt, const std::vector<std::complex<double>>& coefficients) const {
    std::complex<double> cexp = std::exp(1i*cumulative_phase(gt));
    std::vector<std::complex<double>> rho(coefficients);
    // Add back the rotating wave phase to the coherence terms between the
    // low and high states
    for(int k = handler.kmin; k <= handler.kmax; ++k) {
        if(handler.has(1, k, 2, k)) {
            handler.at(rho, 1, k, 2, k) *= cexp;
        }
        if(handler.has(2, k, 1, k)) {
            handler.at(rho, 2, k, 1, k) *= std::conj(cexp);
        }
    }
    return rho;
}

std::vector<std::complex<double>> HMotion::operator()(double gt,
    const std::vector<std::complex<double>>& rho_c) {
    refresh_cache(gt);  // Update cache

    // 1/(i*HBAR) * [H, rho_c] + L(rho_c) from the master equation
    std::vector<std::complex<double>> drho_c(rho_c.size());
#pragma omp parallel for
    for(auto it = handler.idxlist.begin(); it < handler.idxlist.end(); ++it) {
        unsigned nl, nr;
        int kl, kr;
        unsigned idx;
        std::tie(nl, kl, nr, kr, idx) = *it;
        handler.atidx(drho_c, idx) =
            -1i*(haction(rho_c, nl, kl, nr, kr, idx)
                 - std::conj(haction(rho_c, nr, kr, nl, kl)))
            + decayterm(rho_c, nl, kl, nr, kr, idx) * enable_decay;
    }
    return drho_c;
}
