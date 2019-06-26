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
    
    // Set up index maps
    // Use a naive, dense storage for now
    // TODO: change to efficient storage method
    idxmap.reserve(nstates()*nstates());
    idxlist.reserve(idxmap.size());
    for(unsigned nl = 0; nl < nint; ++nl) {
        for(int kl = kmin; kl <= kmax; ++kl) {
            for(unsigned nr = 0; nr < nint; ++nr) {
                for(int kr = kmin; kr <= kmax; ++kr) {
                    idxmap[subidx(nl, kl, nr, kr)] = subidx(nl, kl, nr, kr);
                    idxlist.push_back({nl, kl, nr, kr});
                }
            }
        }
    }
}

unsigned DensMatHandler::stateidx(unsigned n, int k) const {
    // Shift up so the lowest value has index 0
    int k_idx = k - kmin;
    return k_idx + (kmax - kmin + 1)*n;
}
unsigned DensMatHandler::nstates() const {
    return (kmax - kmin + 1)*nint;
}
unsigned DensMatHandler::subidx(unsigned nl, int kl, unsigned nr, int kr) const {
    return stateidx(nr, kr) + nstates()*stateidx(nl, kl);
}

bool DensMatHandler::has(unsigned nl, int kl, unsigned nr, int kr) const {
    return (idxmap.find(subidx(nl, kl, nr, kr)) != idxmap.end());
}

std::complex<double> DensMatHandler::ele(
    const std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    return rho[idxmap.at(subidx(nl, kl, nr, kr))];
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