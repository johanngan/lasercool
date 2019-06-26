#include "HMotion.hpp"
using namespace std::complex_literals;

template<typename T>
inline T sqr(T x) {return x*x;}

HMotion::HMotion(std::string fname):HSwap(fname), SPEED_OF_LIGHT(299792458),
    stationary_decay_prob(0.6), handler(fname) {
    double decay_rate, mass;
    load_params(fname,
        {
            {"spontaneous_decay_rate", &decay_rate},
            {"mass", &mass}
        }
    );

    double k_photon_per_decay = transition_angfreq_per_decay/SPEED_OF_LIGHT;
    recoil_freq_per_decay = HBAR*sqr(k_photon_per_decay)*decay_rate/(2*mass);
}

std::complex<double> HMotion::haction(double gt,
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int kl, unsigned nr, int kr) const {
    std::complex<double> val = 0;

    // Diagonal contribution
    double diag_coeff = recoil_freq_per_decay*sqr(kl);
    switch(nl) {
        case 1: diag_coeff += 0.5*detun_per_decay(gt); break;
        case 2: diag_coeff -= 0.5*detun_per_decay(gt); break;
    }
    val += diag_coeff*handler.ele(rho_c, nl, kl, nr, kr);

    // Off-diagonal contributions
    if(nl > 0) {
        // in rho_c, flip nl: 1 -> 2, 2 -> 1
        unsigned nlflip = !(nl - 1) + 1;
        if(kl - 1 >= handler.kmin) {
            val += 0.5*rabi_softswitch(gt)
                * handler.ele(rho_c, nlflip, kl-1, nr, kr);
        }
        if(kl + 1 <= handler.kmax) {
            val += 0.5*rabi_softswitch(gt)
                * handler.ele(rho_c, nlflip, kl+1, nr, kr);
        }
    }

    return val;
}

std::complex<double> HMotion::decayterm(
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int kl, unsigned nr, int kr) const {
    // On the block diagonal
    if(nl == nr) {
        switch(nl) {
            case 0:
                return (1 - branching_ratio)
                    * handler.ele(rho_c, 2, kl, 2, kr);
            case 1:
            {
                // Approximate anisotropic dipole radiation pattern
                std::complex<double> diprad = stationary_decay_prob
                    * handler.ele(rho_c, 2, kl, 2, kr);
                if(kl-1 >= handler.kmin && kr-1 >= handler.kmin) {
                    diprad += (1-stationary_decay_prob)/2
                        * handler.ele(rho_c, 2, kl-1, 2, kr-1);
                }
                if(kl+1 <= handler.kmax && kr+1 <= handler.kmax) {
                    diprad += (1-stationary_decay_prob)/2
                        * handler.ele(rho_c, 2, kl+1, 2, kr+1);
                }
                return branching_ratio * diprad;
            }
            case 2: // Double decay of coherences within excited state
                return -handler.ele(rho_c, 2, kl, 2, kr);
        }
    } else if(nl == 2 || nr == 2) {
        // Exponential decay of coherences between excited state and lower state
        return -0.5*handler.ele(rho_c, nl, kl, nr, kr);
    }
    return 0;
}

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
    const std::vector<std::complex<double>>& rho_c) const {
    // 1/(i*HBAR) * [H, rho_c] + L(rho_c) from the master equation
    std::vector<std::complex<double>> drho_c(rho_c.size());
    for(const auto& sub: handler.idxlist) {
        unsigned nl = std::get<0>(sub), nr = std::get<2>(sub);
        int kl = std::get<1>(sub), kr = std::get<3>(sub);
        handler.at(drho_c, nl, kl, nr, kr) =
            -1i*haction(gt, rho_c, nl, kl, nr, kr)
            +1i*std::conj(haction(gt, rho_c, nr, kr, nl, kl))
            + decayterm(rho_c, nl, kl, nr, kr) * enable_decay;
    }
    return drho_c;
}