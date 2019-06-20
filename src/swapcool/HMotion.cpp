#include "HMotion.hpp"
using namespace std::complex_literals;

template<typename T>
inline T sqr(T x) {return x*x;}

HMotion::HMotion(std::string fname):HSwap(fname), SPEED_OF_LIGHT(299792458),
    nint(3) {
    double decay_rate, mass, kmax_double;
    load_params(fname,
        {
            {"spontaneous_decay_rate", &decay_rate},
            {"mass", &mass},
            {"max_momentum", &kmax_double}
        }
    );
    kmax = std::abs(static_cast<int>(kmax_double)); // Force to be positive

    double k_photon_per_decay = transition_angfreq_per_decay/SPEED_OF_LIGHT;
    recoil_freq_per_decay = HBAR*sqr(k_photon_per_decay)*decay_rate/(2*mass);
}

unsigned HMotion::stateidx(unsigned n, int k) const {
    // Shift up so the lowest value has index 0
    int k_idx = k + kmax;
    return k_idx + (2*kmax + 1)*n;
}

unsigned HMotion::nstates() const {
    return (2*kmax + 1)*nint;
}

unsigned HMotion::subidx(unsigned nl, int kl, unsigned nr, int kr) const {
    return stateidx(nr, kr) + nstates()*stateidx(nl, kl);
}

unsigned HMotion::nmat() const {
    return nstates()*nstates();
}

std::complex<double> HMotion::totaltr(
    const std::vector<std::complex<double>>& rho_c) const {
    std::complex<double> tr = 0;
    for(unsigned n = 0; n < nint; ++n) {
        for(int k = -kmax; k <= kmax; ++k) {
            tr += rho_c[subidx(n, k, n, k)];
        }
    }
    return tr;
}

std::complex<double> HMotion::partialtr_k(
    const std::vector<std::complex<double>>& rho_c, unsigned n) const {
    std::complex<double> tr = 0;
    for(int k = -kmax; k <= kmax; ++k) {
        tr += rho_c[subidx(n, k, n, k)];
    }
    return tr;
}

std::complex<double> HMotion::partialtr_n(
    const std::vector<std::complex<double>>& rho_c, int k) const {
    std::complex<double> tr = 0;
    for(unsigned n = 0; n < nint; ++n) {
        tr += rho_c[subidx(n, k, n, k)];
    }
    return tr;
}

std::complex<double> HMotion::purity(
    const std::vector<std::complex<double>>& rho_c) const {
    std::complex<double> tr = 0;
    for(unsigned nouter = 0; nouter < nint; ++nouter) {
        for(int kouter = -kmax; kouter <= kmax; ++kouter) {
            for(unsigned ninner = 0; ninner < nint; ++ninner) {
                for(int kinner = -kmax; kinner <= kmax; ++kinner) {
                    tr += rho_c[subidx(nouter, kouter, ninner, kinner)]
                        * rho_c[subidx(ninner, kinner, nouter, kouter)];
                }
            }
        }
    }
    return tr;
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
    val += diag_coeff*rho_c[subidx(nl, kl, nr, kr)];

    // Off-diagonal contributions
    if(nl > 0) {
        // in rho_c, flip nl: 1 -> 2, 2 -> 1
        unsigned nlflip = !(nl - 1) + 1;
        if(kl - 1 >= -kmax) {
            val += 0.5*rabi_softswitch(gt)*rho_c[subidx(nlflip, kl-1, nr, kr)];
        }
        if(kl + 1 <= kmax) {
            val += 0.5*rabi_softswitch(gt)*rho_c[subidx(nlflip, kl+1, nr, kr)];
        }
    }

    return val;
}

std::complex<double> HMotion::decayterm(
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int kl, unsigned nr, int kr) const {
    if(nl == nr && kl == kr) {
        switch(nl) {
            case 0:
                return (1 - branching_ratio)*rho_c[subidx(2, kl, 2, kl)];
            case 1:
            {
                // Approximate anisotropic dipole radiation pattern
                std::complex<double> diprad = 0.6 * rho_c[subidx(2, kl, 2, kl)];
                if(kl - 1 >= -kmax) {
                    diprad += 0.2 * rho_c[subidx(2, kl-1, 2, kl-1)];
                }
                if(kl + 1 <= kmax) {
                    diprad += 0.2 * rho_c[subidx(2, kl+1, 2, kl+1)];
                }
                return branching_ratio * diprad;
            }
            case 2:
                return -rho_c[subidx(2, kl, 2, kl)];
        }
    } else if(nl == 2 || nr == 2) {
        // Exponential decay of coherences with excited state
        return -0.5*rho_c[subidx(nl, kl, nr, kr)];
    }
    return 0;
}

std::vector<std::complex<double>> HMotion::density_matrix(
    double gt, const std::vector<std::complex<double>>& coefficients) const {
    std::complex<double> cexp = std::exp(1i*cumulative_phase(gt));
    std::vector<std::complex<double>> rho(coefficients);
    // Add back the rotating wave phase to the coherence terms between the
    // low and high states
    for(int k = -kmax; k <= kmax; ++k) {
        rho[subidx(1, k, 2, k)] *= cexp;
        rho[subidx(2, k, 1, k)] *= std::conj(cexp);
    }
    return rho;
}

std::vector<std::complex<double>> HMotion::operator()(double gt,
    const std::vector<std::complex<double>>& rho_c) const {
    // 1/(i*HBAR) * [H, rho_c] + L(rho_c) from the master equation
    std::vector<std::complex<double>> drho_c(rho_c.size());
    for(unsigned nl = 0; nl < nint; ++nl) {
        for(int kl = -kmax; kl <= kmax; ++kl) {
            for(unsigned nr = 0; nr < nint; ++nr) {
                for(int kr = -kmax; kr <= kmax; ++kr) {
                    drho_c[subidx(nl, kl, nr, kr)] =
                        -1i*haction(gt, rho_c, nl, kl, nr, kr)
                        +1i*std::conj(haction(gt, rho_c, nr, kr, nl, kl))
                        + decayterm(rho_c, nl, kl, nr, kr) * enable_decay;
                }
            }
        }
    }
    return drho_c;
}