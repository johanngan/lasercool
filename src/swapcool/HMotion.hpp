#ifndef HMOTION_HPP_
#define HMOTION_HPP_

#include "HSwap.hpp"

// Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency, under the rotating wave approximation,
// including interaction with the laser and also motional states
struct HMotion : public HSwap {
    const double SPEED_OF_LIGHT;
    unsigned nint;  // number of internal states
    int kmax;   // maximum k value (assumed positive)
    double recoil_freq_per_decay;

    HMotion(std::string);
    // Convert state subscripts to linear indexes, enumerated as |n, k>
    unsigned stateidx(unsigned, int) const;
    // Number of entries in a state
    unsigned nstates() const;
    // Convert state subscripts to linear indexes in the density matrix,
    // enumerated as |n-left, k-left><n-right, k-right|
    unsigned subidx(unsigned, int, unsigned, int) const;
    // Number of entries in the matrix
    unsigned nmat() const;

    // Total trace
    std::complex<double> totaltr(const std::vector<std::complex<double>>&) const;

    // Partial trace over k for a fixed n
    std::complex<double> partialtr_k(const std::vector<std::complex<double>>&,
        unsigned) const;

    // Partial trace over n for a fixed k
    std::complex<double> partialtr_n(const std::vector<std::complex<double>>&,
        int) const;
    
    // Trace of rho^2
    std::complex<double> purity(const std::vector<std::complex<double>>&) const;

    // The action of the Hamiltonian on the density matrix, returns a single
    // component of H*rho, enumerated as in subidx()
    std::complex<double> haction(double,
        const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int) const;
    
    // The spontaneous decay part of the derivative (Lindblad superoperator)
    std::complex<double> decayterm(const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int) const;

    // Transforms the coefficients solved for in the rotating wave
    // approximation back to the actual density matrix values;
    // i.e. put the oscillation back in.
    std::vector<std::complex<double>> density_matrix(
        double, const std::vector<std::complex<double>>&) const override;

    // Derivative operator to be passed to the timestepper
    std::vector<std::complex<double>> operator()(double,
        const std::vector<std::complex<double>>&) const override;
};

#endif