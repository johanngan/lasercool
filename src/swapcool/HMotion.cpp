#include "HMotion.hpp"
using namespace std::complex_literals;

template<typename T>
inline T sqr(T x) {return x*x;}

// Analytical form of the dipole integral between theta1 and theta2
// Integral from theta1 to theta2 of 2/pi*sin^2(theta)d(theta)
double dipole_integral(double theta1, double theta2) {
    return (theta2 - theta1 - 0.5*(sin(2*theta2) - sin(2*theta1))) / M_PI;
}

HMotion::HMotion(std::string fname):HSwap(fname) {
    double mass, init_temp, ksigmas, kmin_double, kmax_double, ksubdivs_double;
    load_params(fname,
        {
            {"mass", &mass},
            {"initial_temperature", &init_temp},
            {"momentum_stddevs", &ksigmas},
            {"min_momentum", &kmin_double},
            {"max_momentum", &kmax_double},
            {"momentum_subdivisions", &ksubdivs_double}
        }
    );
    // Default to 1 subdivision (i.e. only tracking integer values)
    if(std::isnan(ksubdivs_double)) {
        ksubdivs_double = 1;
    }

    double k_photon_per_decay = transition_angfreq_per_decay
        /fundamental_constants::SPEED_OF_LIGHT;
    recoil_freq_per_decay = fundamental_constants::HBAR
        *sqr(k_photon_per_decay)*decay_rate/(2*mass);

    int kmin, kmax;
    // Manually set k range
    if(!std::isnan(kmax_double)) {
        if(std::isnan(kmin_double)) {
            kmin_double = -kmax_double;
        }
        kmin = static_cast<int>(kmin_double);
        kmax = static_cast<int>(kmax_double);
    } else {
        if(std::isnan(ksigmas)) {
            // Default to 3 standard deviations away
            ksigmas = 3;
        }
        // Calculate k range from standard deviation
        kmax = static_cast<int>(std::round(ksigmas*
            sqrt(fundamental_constants::K_BOLTZMANN*init_temp
                / (2*fundamental_constants::HBAR*recoil_freq_per_decay*decay_rate))
        ));
        kmin = -kmax;
    }

    // Subdivide the k states
    handler = DensMatHandler(kmin, kmax, static_cast<int>(ksubdivs_double));

    // Compute the dipole radiation probabilities for each sector
    diprad_probs.reserve(handler.ksubdivs+1);
    for(int dkidx = 0; dkidx <= handler.ksubdivs; ++dkidx) {
        // The actual k value is the k index divided by the number of
        // subdivisions, m.
        // Sector for each k value has boundaries of theta where
        // cos(theta) = (k +- 1/2)/m
        // halfway between the sector's k value and those of the sectors
        // immediately above and below. The boundaries are capped at
        // k/m = 0, 1, so as not to overcount k/m = 1(forward) and to
        // half-count k/m = 0 (transverse), to compensate for double counting
        // when implementing the population decay term.
        double upper = std::min(1., (dkidx + 0.5) / handler.ksubdivs);
        double lower = std::max(0., (dkidx - 0.5) / handler.ksubdivs);
        diprad_probs.push_back(dipole_integral(
            acos(upper), acos(lower)));
    }
}

std::complex<double> HMotion::haction(
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int kl, unsigned nr, int kr, int idx) const {
    // Read from cache
    double cachehalfdetun = cache[halfdetun], cachehalfrabi = cache[halfrabi];

    std::complex<double> val = 0;

    // Diagonal contribution
    double diag_coeff = recoil_freq_per_decay*sqr(handler.kval(kl));
    switch(nl) {
        case 1: diag_coeff += cachehalfdetun; break;
        case 2: diag_coeff -= cachehalfdetun; break;
    }
    if(idx != -1) {
        // Use precomputed index
        val += diag_coeff*handler.atidx(rho_c, idx);
    } else {
        val += diag_coeff*handler.ele(rho_c, nl, kl, nr, kr);
    }

    // Off-diagonal contributions
    if(nl > 0) {
        // in rho_c, flip nl: 1 -> 2, 2 -> 1
        unsigned nlflip = !(nl - 1) + 1;
        // Kick the index up and down by ksubdivs, so that the actual k value
        // gets a full kick of dk = 1
        if(kl - handler.ksubdivs >= handler.kmin) {
            val += cachehalfrabi*handler.ele(rho_c,
                nlflip, kl-handler.ksubdivs, nr, kr);
        }
        if(kl + handler.ksubdivs <= handler.kmax) {
            val += cachehalfrabi*handler.ele(rho_c,
                nlflip, kl+handler.ksubdivs, nr, kr);
        }
    }

    return val;
}

std::complex<double> HMotion::decayterm(
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int kl, unsigned nr, int kr, unsigned idx) const {
    // On the block diagonal
    if(nl == nr) {
        switch(nl) {
            case 0:
                return (1 - branching_ratio)
                    * handler.ele(rho_c, 2, kl, 2, kr);
            case 1:
            {
                // Approximate anisotropic dipole radiation pattern
                std::complex<double> diprad = 0;
                // Double counts dk = 0, but the diprad_probs entry is halved, so
                // it evens out
                for(int dkidx = 0; dkidx < static_cast<int>(diprad_probs.size());
                    ++dkidx) {
                    // Decay from below
                    if(kl - dkidx >= handler.kmin && kr - dkidx >= handler.kmin) {
                        diprad += diprad_probs[dkidx] * handler.ele(rho_c,
                            2, kl - dkidx, 2, kr - dkidx);
                    }
                    // Decay from above
                    if(kl + dkidx <= handler.kmax && kr + dkidx <= handler.kmax) {
                        diprad += diprad_probs[dkidx] * handler.ele(rho_c,
                            2, kl + dkidx, 2, kr + dkidx);
                    }
                }
                return branching_ratio * diprad;
            }
            case 2: // Double decay of coherences within excited state
                // idx is guaranteed to be stored by design
                return -handler.atidx(rho_c, idx);
        }
    } else if(nl == 2 || nr == 2) {
        // Exponential decay of coherences between excited state and lower state
        // idx is guaranteed to be stored by design
        return -0.5*handler.atidx(rho_c, idx);
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

void HMotion::initialize_cycle(std::vector<std::complex<double>>& rho) const {
    for(int kl = handler.kmin; kl <= handler.kmax; ++kl) {
        for(int kr = handler.kmin; kr <= handler.kmax; ++kr) {
            // Excited state population and intra-excited-state coherences
            // distribute between the lower energy states
            if(handler.has(0, kl, 0, kr)) {
                handler.at(rho, 0, kl, 0, kr) += (1 - branching_ratio)
                    * handler.ele(rho, 2, kl, 2, kr);
            }

            // Double counts dk = 0, but the diprad_probs entry is halved, so
            // it evens out
            for(int dkidx = 0; dkidx < static_cast<int>(diprad_probs.size());
                ++dkidx) {
                // Decay up
                if(handler.has(1, kl + dkidx, 1, kr + dkidx)) {
                    if(kl + dkidx <= handler.kmax && kr + dkidx <= handler.kmax) {
                        handler.at(rho, 1, kl + dkidx, 1, kr + dkidx) +=
                            diprad_probs[dkidx] * handler.ele(rho, 2, kl, 2, kr);
                    }
                }
                // Decay down
                if(handler.has(1, kl - dkidx, 1, kr - dkidx)) {
                    if(kl - dkidx >= handler.kmin && kr - dkidx >= handler.kmin) {
                        handler.at(rho, 1, kl - dkidx, 1, kr - dkidx) +=
                            diprad_probs[dkidx] * handler.ele(rho, 2, kl, 2, kr);
                    }
                }
            }

            // Excited state and excited-state coherences decay to 0
            if(handler.has(2, kl, 2, kr)) {
                handler.at(rho, 2, kl, 2, kr) = 0;
            }
            if(handler.has(1, kl, 2, kr)) {
                handler.at(rho, 1, kl, 2, kr) = 0;
            }
            if(handler.has(2, kl, 1, kr)) {
                handler.at(rho, 2, kl, 1, kr) = 0;
            }
        }
    }
    return;
}