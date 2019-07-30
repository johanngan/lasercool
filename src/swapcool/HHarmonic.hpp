#ifndef HHARMONIC_HPP_
#define HHARMONIC_HPP_

#include <exception>
#include <unordered_map>
#include <algorithm>
#include "HMotion.hpp"
#include "lasercool/fundconst.hpp"

// SWAP Hamiltonian for simple harmonic oscillator states
struct HHarmonic : public HMotion {
    // Oscillator omega per decay rate
    double osc_angfreq_per_decay;
    int coupling_limit;
    // Couplings between oscillator states due to the driving lasers
    // Outer (first) key corresponds to destination state, inner (second) key
    // corresponds to source state, value corresponds to overlap amplitude from
    // cos(kx)
    std::unordered_map<int, std::unordered_map<int, double>> drive_couplings;
    // Couplings between oscillator states due to spontaneous decay
    // Outer (first) key corresponds to destination state, inner (second) key
    // corresponds to source state, value corresponds to overlap amplitudes
    // from e^(-ikx)
    std::unordered_map<int, std::unordered_map<int, std::complex<double>>>
        decay_couplings;
    // Couplings of density matrix indexes for other states state decaying into one
    std::unordered_map<int, std::unordered_map<int, std::vector<
        std::tuple<int, int, std::complex<double>>>>> decay_inverse_couplings;

    HHarmonic(std::string);

    // Returns a map of the coefficients of each |h'> when (a + a^dag)^power is
    // applied to |h>,
    // Computes coefficient recursively from lower power results, and stores them
    // in a cache
    static std::unordered_map<int, double> x_pwr_p_operator(int, int,
        std::unordered_map<int, std::unordered_map<int, double>>&);
    // Action on |h> of an operator that's some exponential-like power
    // series in (prefactor*(a+a^dag))^p, i.e. terms are of the form +-x^p/(p!)
    static std::unordered_map<int, std::complex<double>> x_power_series_operator(
        int, std::complex<double>, bool, int, int, int,
        std::unordered_map<int, std::unordered_map<int, double>>&);
    // Returns a map of the coefficients of each |h'> when cos(kx) is applied to |h>
    static std::unordered_map<int, double> cosx_operator(int, double, int,
        std::unordered_map<int, std::unordered_map<int, double>>&);
    // Returns a map of the coefficients of each |h'> when e^(ikx) is applied to |h>
    static std::unordered_map<int, std::complex<double>> iexpx_operator(int,
        double, int, std::unordered_map<int, std::unordered_map<int, double>>&);

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