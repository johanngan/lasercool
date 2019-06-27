#include "DensMatHandler.hpp"

DensMatHandler::DensMatHandler(std::string fname):nint(3) {
    double kmax_double, kmin_double;
    load_params(fname,
        {
            {"max_momentum", &kmax_double},
            {"min_momentum", &kmin_double}
        }
    );
    // Default value
    if(std::isnan(kmin_double)) {
        kmin_double = -kmax_double;
    }
    kmax = static_cast<int>(kmax_double);
    kmin = static_cast<int>(kmin_double);
    if(kmin > kmax) {
        throw std::runtime_error("Min momentum greater than max momentum.");
    }
    // Calculate state numbers
    kstates = kmax - kmin + 1;
    
    // Set up index maps
    // Store only the upper triangle, and also exclude the coherences with
    // the ground state
    idxmap.reserve(3*kstates*(kstates+1)/2 + kstates*kstates);
    idxlist.reserve(idxmap.size());
    // On the upper triangles of the block diagonal
    for(unsigned n = 0; n < nint; ++n) {
        for(int kl = kmin; kl <= kmax; ++kl) {
            for(int kr = kl; kr <= kmax; ++kr) {
                idxmap[subidx(n, kl, n, kr)] = idxmap.size();
                idxlist.push_back({n, kl, n, kr});
            }
        }
    }
    // The upper-triangular coherences between the high and low energy states
    for(int kl = kmin; kl <= kmax; ++kl) {
        for(int kr = kmin; kr <= kmax; ++kr) {
            idxmap[subidx(1, kl, 2, kr)] = idxmap.size();
            idxlist.push_back({1, kl, 2, kr});
        }
    }
}

inline unsigned DensMatHandler::subidx(
    unsigned nl, int kl, unsigned nr, int kr) const {
    return kr-kmin + kstates*(nr + nint*(kl-kmin + kstates*nl));
}

bool DensMatHandler::has(unsigned nl, int kl, unsigned nr, int kr) const {
    return (idxmap.find(subidx(nl, kl, nr, kr)) != idxmap.end());
}

std::complex<double> DensMatHandler::ele(
    const std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    // Directly stored
    if(has(nl, kl, nr, kr)) {
        return rho[idxmap.at(subidx(nl, kl, nr, kr))];
    } else if(has(nr, kr, nl, kl)) { // Try the transpose
        // Density matrix must be Hermitian
        return std::conj(rho[idxmap.at(subidx(nr, kr, nl, kl))]);
    }
    // If nothing is found, must be a 0 entry (coherence with ground state)
    return 0;
}
std::complex<double>& DensMatHandler::at(
    std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    return rho[idxmap.at(subidx(nl, kl, nr, kr))];
}
const std::complex<double>& DensMatHandler::at(
    const std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    return rho[idxmap.at(subidx(nl, kl, nr, kr))];
}

std::complex<double> DensMatHandler::totaltr(
    const std::vector<std::complex<double>>& rho_c) const {
    std::complex<double> tr = 0;
    for(unsigned n = 0; n < nint; ++n) {
        for(int k = kmin; k <= kmax; ++k) {
            tr += ele(rho_c, n, k, n, k);
        }
    }
    return tr;
}

std::complex<double> DensMatHandler::partialtr_k(
    const std::vector<std::complex<double>>& rho_c, unsigned n) const {
    std::complex<double> tr = 0;
    for(int k = kmin; k <= kmax; ++k) {
        tr += ele(rho_c, n, k, n, k);
    }
    return tr;
}

std::complex<double> DensMatHandler::partialtr_n(
    const std::vector<std::complex<double>>& rho_c, int k) const {
    std::complex<double> tr = 0;
    for(unsigned n = 0; n < nint; ++n) {
        tr += ele(rho_c, n, k, n, k);
    }
    return tr;
}

std::complex<double> DensMatHandler::purity(
    const std::vector<std::complex<double>>& rho_c) const {
    std::complex<double> tr = 0;
    for(unsigned nouter = 0; nouter < nint; ++nouter) {
        for(int kouter = kmin; kouter <= kmax; ++kouter) {
            for(unsigned ninner = 0; ninner < nint; ++ninner) {
                for(int kinner = kmin; kinner <= kmax; ++kinner) {
                    tr += ele(rho_c, nouter, kouter, ninner, kinner)
                        * ele(rho_c, ninner, kinner, nouter, kouter);
                }
            }
        }
    }
    return tr;
}