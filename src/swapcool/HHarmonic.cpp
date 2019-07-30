#include "HHarmonic.hpp"
using namespace std::complex_literals;

template<typename T>
inline T sqr(T x) {return x*x;}

// Add population to a couplings map given action of an operator
template<typename T>
inline void add_couplings(
    std::unordered_map<int, std::unordered_map<int, T>>& couplings_map,
    std::unordered_map<int, T> operator_map, int h, int hmin, int hmax,
    int coupling_limit) {
    for(const auto op_action: operator_map) {
        int state = op_action.first;
        T coeff = op_action.second;
        // Only add if the source state is within the basis
        if(state >= hmin && state <= hmax &&
            std::abs(state - h) <= coupling_limit) {
            couplings_map[state][h] += coeff;
        }
    }
}

HHarmonic::HHarmonic(std::string fname):HMotion(fname) {
    double mass, init_temp, osc_angfreq, hdecays, hmin_double, hmax_double,
        max_order_double, coupling_limit_double;
    load_params(fname,
        {
            {"mass", &mass},
            {"initial_temperature", &init_temp},
            {"oscillator_angular_frequency", &osc_angfreq},
            {"energy_decay_constants", &hdecays},
            {"min_energy_state", &hmin_double},
            {"max_energy_state", &hmax_double},
            {"series_truncation_order", &max_order_double},
            {"coupling_limit", &coupling_limit_double}
        }
    );
    osc_angfreq_per_decay = osc_angfreq / decay_rate;
    double k_photon = transition_angfreq_per_decay*decay_rate
        /fundamental_constants::SPEED_OF_LIGHT;
    // Effective oscillator amplitude in the classical analogy
    double osc_amplitude = sqrt(fundamental_constants::HBAR/(mass*osc_angfreq));
    // Coefficient for the position operator in front of the ladder operators
    // for the kx operator
    double ladder_coeff = k_photon*osc_amplitude/sqrt(2);

    // Manually set oscillator state range
    // hmin defaults to 0
    int hmin = 0, hmax;
    if(!std::isnan(hmin_double)) {
        hmin = static_cast<int>(hmin_double);
        if(hmin < 0) {
            throw std::invalid_argument("Cannot have a negative oscillator state.");
        }
    }

    if(!std::isnan(hmax_double)) {
        hmax = static_cast<int>(hmax_double);
    } else {
        // Default to 5 decay constants away from h = 0
        if(std::isnan(hdecays)) {
            hdecays = 5;
        }
        // Calculate h states from decay constants
        hmax = static_cast<int>(std::round(hdecays*
            fundamental_constants::K_BOLTZMANN*init_temp
            / (fundamental_constants::HBAR*osc_angfreq)));
    }
    handler = DensMatHandler(hmin, hmax);
    int hstates = hmax - hmin + 1;

    // Default to a truncation order of 5*sqrt(hmax) terms (but at least 20),
    // since the amplitudes asymptotically go like h^(dh/2) / (dh!)
    // (up to and including x^(max_order) term)
    int max_order = std::max(20, static_cast<int>(std::ceil(5*sqrt(hmax))));
    if(!std::isnan(max_order_double)) {
        max_order = static_cast<int>(max_order_double)-5;
    }
    if(!std::isnan(coupling_limit_double)) {
        // Coupling can't go past the max order
        coupling_limit = std::min(max_order, static_cast<int>(coupling_limit_double));
    } else {
        coupling_limit = max_order;
    }

    // Precompute the couplings between oscillator states due to the laser and
    // spontaneous decay
    drive_couplings.reserve(hstates);
    decay_couplings.reserve(hstates);
    for(int h = hmin; h <= hmax; ++h) {
        std::unordered_map<int, std::unordered_map<int, double>> coeffs_cache;
        // States will couple up and down by max_order, but can't go below 0
        coeffs_cache.reserve(hmax+max_order - std::max(0, (hmin-max_order)) + 1);

        // Compute decay couplings from the exponential (-)
        add_couplings(decay_couplings,
            iexpx_operator(h, -ladder_coeff, max_order, coeffs_cache),
            h, hmin, hmax, coupling_limit);

        // Compute drive couplings from cosine
        add_couplings(drive_couplings,
            cosx_operator(h, ladder_coeff, max_order, coeffs_cache),
            h, hmin, hmax, coupling_limit);
    }

    // Precompute composite decay couplings for the density matrix
    decay_inverse_couplings.reserve(hstates);
    for(int hl = handler.kmin; hl <= handler.kmax; ++hl) {
        decay_inverse_couplings[hl].reserve(hstates);
        for(int hr = handler.kmin; hr <= handler.kmax; ++hr) {
            int hlothermin = std::max(handler.kmin, hl-coupling_limit),
                hlothermax = std::min(handler.kmax, hl+coupling_limit),
                hrothermin = std::max(handler.kmin, hr-coupling_limit),
                hrothermax = std::min(handler.kmax, hr+coupling_limit);
            int nreserve = (hlothermax-hlothermin+1)*(hrothermax-hrothermin+1);
            decay_inverse_couplings[hl][hr].reserve(nreserve);
            for(int hlother = hlothermin; hlother <= hlothermax; ++hlother) {
                for(int hrother = hrothermin; hrother <= hrothermax; ++hrother) {
                    // Contribution from the left component decaying from hlother
                    // to hl, and the right component decaying from hrother to hr;
                    // contributions from e^(-ikx) (right kick)
                    // and e^(ikx) (left kick)
                    std::complex<double> inv_couple = 
                        // A kick only happens with some probability
                        (1 - stationary_decay_prob) / 2 * (
                            // From right kick decays
                            decay_couplings[hl][hlother]
                                *std::conj(decay_couplings[hr][hrother]) +
                            // From left kick decays
                            std::conj(decay_couplings[hlother][hl])
                                *decay_couplings[hrother][hr]
                        );
                    decay_inverse_couplings[hl][hr].push_back(
                        std::make_tuple(hlother, hrother, inv_couple));
                }
            }
        }
    }
}

std::unordered_map<int, double> HHarmonic::x_pwr_p_operator(int h, int power, 
    std::unordered_map<int, std::unordered_map<int, double>>& coeffs_cache) {
    // Check cache
    if(coeffs_cache.find(power) != coeffs_cache.end()) {
        return coeffs_cache[power];
    }
    // Base case; identity operator
    if(power == 0) {
        coeffs_cache[power][h] = 1;
        return coeffs_cache[power];
    }
    // Recursively fill the cache up to the previous power
    x_pwr_p_operator(h, power - 1, coeffs_cache);
    // Go through each final state of the previous power to build the next
    std::unordered_map<int, double> coeffs;
    for(const auto& prev: coeffs_cache[power - 1]) {
        int prev_state = prev.first;
        double prev_coeff = prev.second;
        // Apply the raising operator
        coeffs[prev_state + 1] += prev_coeff * sqrt(prev_state + 1);
        // Apply the lowering operator
        if(prev_state > 0) {
            coeffs[prev_state - 1] += prev_coeff * sqrt(prev_state);
        }
    }
    coeffs_cache[power] = coeffs;
    return coeffs;
}

std::unordered_map<int, std::complex<double>> HHarmonic::x_power_series_operator(
    int h, std::complex<double> prefactor, bool alternating,
    int min_order, int max_order, int order_stride,
    std::unordered_map<int, std::unordered_map<int, double>>& coeffs_cache) {
    std::unordered_map<int, std::complex<double>> coeffs;
    std::complex<double> taylor_coeff = 1;    // Taylor coefficient
    for(int order = 0; order <= max_order; ++order) {
        // Only add to the coefficients every stride, after min_order is reached
        if(order >= min_order && (order - min_order) % order_stride == 0) {
            for(const auto& p_coeffs: x_pwr_p_operator(
                h, order, coeffs_cache)) {
                coeffs[p_coeffs.first] += taylor_coeff * p_coeffs.second;
            }
            if(alternating) {
                taylor_coeff *= -1; // Alternate sign
            }
        }
        // Build up the Taylor coefficient every order, even if not included
        // in the series
        // numerator does the exponent, denominator does the factorial
        taylor_coeff *= prefactor / std::complex<double>(order + 1);
    }
    return coeffs;
}

// cos(xprefactor * (a+a^dag))
std::unordered_map<int, double> HHarmonic::cosx_operator(
    int h, double xprefactor, int max_order,
    std::unordered_map<int, std::unordered_map<int, double>>& coeffs_cache) {
    auto coeffs_complex = x_power_series_operator(h, xprefactor, true, 0,
        max_order, 2, coeffs_cache);
    std::unordered_map<int, double> coeffs;
    coeffs.reserve(coeffs_complex.size());
    for(const auto& coeff: coeffs_complex) {
        coeffs[coeff.first] = std::real(coeff.second);
    }
    return coeffs;
}

// e^(i*xprefactor*(a+a^dag)), magnitude squared of the amplitudes
std::unordered_map<int, std::complex<double>> HHarmonic::iexpx_operator(
    int h, double xprefactor, int max_order,
    std::unordered_map<int, std::unordered_map<int, double>>& coeffs_cache) {
    return x_power_series_operator(h, 1i*xprefactor, false, 0, max_order, 1,
        coeffs_cache);
}

std::complex<double> HHarmonic::haction(
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int hl, unsigned nr, int hr, int idx) const {
    // Read from cache
    double cachehalfdetun = cache[halfdetun], cachehalfrabi = cache[halfrabi];

    std::complex<double> val = 0;

    // Diagonal contribution
    double diag_coeff = (hl + 0.5)*osc_angfreq_per_decay;
    switch(nl) {
        case 1: diag_coeff += cachehalfdetun; break;
        case 2: diag_coeff -= cachehalfdetun; break;
    }
    if(idx != -1) {
        // Use precomputed index
        val += diag_coeff*handler.atidx(rho_c, idx);
    } else {
        val += diag_coeff*handler.ele(rho_c, nl, hl, nr, hr);
    }

    // Off-diagonal contributions
    if(nl > 0) {
        // in rho_c, flip nl: 1 -> 2, 2 -> 1
        unsigned nlflip = !(nl - 1) + 1;

        // Sum over all contributions to this state coming from all other states,
        // weighted by their overlap coefficients
        for(const auto& coupling: drive_couplings.at(hl)) {
            val += cachehalfrabi * coupling.second * handler.ele(
                rho_c, nlflip, coupling.first, nr, hr);
        }
    }

    return val;
}

std::complex<double> HHarmonic::decayterm(
    const std::vector<std::complex<double>>& rho_c,
    unsigned nl, int hl, unsigned nr, int hr, unsigned idx) const {
    // On the block diagonal
    if(nl == nr) {
        switch(nl) {
            case 0:
                return (1 - branching_ratio)
                    * handler.ele(rho_c, 2, hl, 2, hr);
            case 1:
            {
                // No kick contribution
                std::complex<double> diprad = stationary_decay_prob *
                    handler.ele(rho_c, 2, hl, 2, hr);
                // Sum over all the matrix elements that have decay contributions
                // into this matrix element
                for(const auto& source: decay_inverse_couplings.at(hl).at(hr)) {
                    int hl0, hr0;
                    std::complex<double> coupling;
                    std::tie(hl0, hr0, coupling) = source;
                    diprad += coupling * handler.ele(rho_c, 2, hl0, 2, hr0);
                }
                // Decay to state 1 only happens with some branching ratio
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

void HHarmonic::initialize_cycle(std::vector<std::complex<double>>& rho) const {
    // Only run decays if they're enabled
    if(!enable_decay) return;

    // Accumulate into state 1
    for(int hl = handler.kmin; hl <= handler.kmax; ++hl) {
        for(int hr = handler.kmin; hr <= handler.kmax; ++hr) {
            // Excited state population and intra-excited-state coherences
            // distribute between the lower energy states
            if(handler.has(0, hl, 0, hr)) {
                handler.at(rho, 0, hl, 0, hr) += (1 - branching_ratio)
                    * handler.ele(rho, 2, hl, 2, hr);
            }
            // Stationary decay
            if(handler.has(1, hl, 1, hr)) {
                handler.at(rho, 1, hl, 1, hr) +=
                    stationary_decay_prob*branching_ratio
                    * handler.ele(rho, 2, hl, 2, hr);
            }
            // Decay from kicks
            if(handler.has(1, hl, 1, hr)) {
                // Sum over all the matrix elements that have decay contributions
                // into this matrix element
                for(const auto& source: decay_inverse_couplings.at(hl).at(hr)) {
                    int hl0, hr0;
                    std::complex<double> coupling;
                    std::tie(hl0, hr0, coupling) = source;
                    handler.at(rho, 1, hl, 1, hr) +=
                        branching_ratio*coupling*handler.ele(rho, 2, hl0, 2, hr0);
                }
            }
        }
    }

    // Decay away state 2
    for(int hl = handler.kmin; hl <= handler.kmax; ++hl) {
        for(int hr = handler.kmin; hr <= handler.kmax; ++hr) {
            // Excited state and excited-state coherences decay to 0
            if(handler.has(2, hl, 2, hr)) {
                handler.at(rho, 2, hl, 2, hr) = 0;
            }
            if(handler.has(1, hl, 2, hr)) {
                handler.at(rho, 1, hl, 2, hr) = 0;
            }
            if(handler.has(2, hl, 1, hr)) {
                handler.at(rho, 2, hl, 1, hr) = 0;
            }
        }
    }
}