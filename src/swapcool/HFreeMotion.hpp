#ifndef HFREEMOTION_HPP_
#define HFREEMOTION_HPP_

#include "HMotion.hpp"
#include "lasercool/fundconst.hpp"

// SWAP Hamiltonian with free-particle momentum states
struct HFreeMotion : public HMotion {
    double recoil_freq_per_decay;

    HFreeMotion(std::string);

    // The action of the Hamiltonian on the density matrix, returns a single
    // component of H*rho
    // Optionally provide a precomputed stored index for speed
    // -1 means no index is provided
    std::complex<double> haction(const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int, int idx=-1) const override;
    
    // The spontaneous decay part of the derivative (Lindblad superoperator)
    std::complex<double> decayterm(const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int, unsigned) const override;

    // Modify the density matrix in preparation for a new cycle
    void initialize_cycle(std::vector<std::complex<double>>&) const override;
};

#endif