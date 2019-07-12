#include "DensMatHandler.hpp"

DensMatHandler::DensMatHandler(int kvalmin, int kvalmax, int ksubdivs):
    nint(3), kmin(ksubdivs*kvalmin), kmax(ksubdivs*kvalmax), ksubdivs(ksubdivs) {
    if(kmin > kmax) {
        throw std::invalid_argument("Min momentum greater than max momentum.");
    }
    // Calculate state numbers/increments
    kstates = kmax - kmin + 1;
    krinc = 1;
    nrinc = kstates*krinc;
    klinc = nint*nrinc;
    nlinc = kstates*klinc;
    ninc = nlinc + nrinc;
    kinc = klinc + krinc;
    
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
                idxlist.push_back({n, kl, n, kr, subidx(n, kl, n, kr)});
            }
        }
    }
    // The upper-triangular coherences between the high and low energy states
    for(int kl = kmin; kl <= kmax; ++kl) {
        for(int kr = kmin; kr <= kmax; ++kr) {
            idxmap[subidx(1, kl, 2, kr)] = idxmap.size();
            idxlist.push_back({1, kl, 2, kr, subidx(1, kl, 2, kr)});
        }
    }
}

double DensMatHandler::kval(int kidx) const {
    return static_cast<double>(kidx) / ksubdivs;
}

int DensMatHandler::closest_kidx(double kval) const {
    return std::round(kval * ksubdivs);
}

unsigned DensMatHandler::subidx(
    unsigned nl, int kl, unsigned nr, int kr) const {
    return kr-kmin + kstates*(nr + nint*(kl-kmin + kstates*nl));
}

bool DensMatHandler::has(unsigned nl, int kl, unsigned nr, int kr) const {
    return hasidx(subidx(nl, kl, nr, kr));
}
bool DensMatHandler::hasidx(unsigned idx) const {
    return (idxmap.find(idx) != idxmap.end());
}

std::complex<double> DensMatHandler::ele(
    const std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    return eleidx(rho, nl, kl, nr, kr, subidx(nl, kl, nr, kr));
}
std::complex<double> DensMatHandler::eleidx(
    const std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr, unsigned idx) const {
    // Directly stored
    if(hasidx(idx)) {
        return atidx(rho, idx);
    } else {
        // Try the transpose
        unsigned idxtrans = idxtranspose(nl, kl, nr, kr, idx);
        if(hasidx(idxtrans)) {
            // Density matrix must be Hermitian
            return std::conj(atidx(rho, idxtrans));
        }
    }
    // If nothing is found, must be a 0 entry (coherence with ground state)
    return 0;
}
std::complex<double>& DensMatHandler::at(
    std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    return atidx(rho, subidx(nl, kl, nr, kr));
}
const std::complex<double>& DensMatHandler::at(
    const std::vector<std::complex<double>>& rho,
    unsigned nl, int kl, unsigned nr, int kr) const {
    return atidx(rho, subidx(nl, kl, nr, kr));
}
std::complex<double>& DensMatHandler::atidx(
    std::vector<std::complex<double>>& rho, unsigned idx) const {
    return rho[idxmap.at(idx)];
}
const std::complex<double>& DensMatHandler::atidx(
    const std::vector<std::complex<double>>& rho, unsigned idx) const {
    return rho[idxmap.at(idx)];
}

unsigned DensMatHandler::swapsub(int sub1, int sub2, int inc1, int inc2,
    unsigned idx) const {
    return idx + (inc1 - inc2)*(sub2 - sub1);
}
unsigned DensMatHandler::idxtranspose(unsigned nl, int kl, unsigned nr, int kr,
    unsigned idx) const {
    return swapsub(kl, kr, klinc, krinc,
        swapsub(nl, nr, nlinc, nrinc, idx));
}

std::complex<double> DensMatHandler::totaltr(
    const std::vector<std::complex<double>>& rho_c) const {
    std::complex<double> tr = 0;
    for(unsigned n = 0; n < nint; ++n) {
        for(int kidx = kmin; kidx <= kmax; ++kidx) {
            tr += ele(rho_c, n, kidx, n, kidx);
        }
    }
    return tr;
}

std::complex<double> DensMatHandler::partialtr_k(
    const std::vector<std::complex<double>>& rho_c, unsigned n) const {
    std::complex<double> tr = 0;
    for(int kidx = kmin; kidx <= kmax; ++kidx) {
        tr += ele(rho_c, n, kidx, n, kidx);
    }
    return tr;
}

std::complex<double> DensMatHandler::partialtr_n(
    const std::vector<std::complex<double>>& rho_c, int kidx) const {
    std::complex<double> tr = 0;
    for(unsigned n = 0; n < nint; ++n) {
        tr += ele(rho_c, n, kidx, n, kidx);
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