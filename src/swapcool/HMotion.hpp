#ifndef HMOTION_HPP_
#define HMOTION_HPP_

#include "HSwap.hpp"
#include "DensMatHandler.hpp"

// Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency, under the rotating wave approximation,
// including interaction with the laser and also motional states
struct HMotion : public HSwap {
    const double SPEED_OF_LIGHT;
    // probability to decay from excited state without changing momentum
    double stationary_decay_prob;
    double recoil_freq_per_decay;
    DensMatHandler handler;

    HMotion(std::string);

    // The action of the Hamiltonian on the density matrix, returns a single
    // component of H*rho
    std::complex<double> haction(double,
        const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int);
    
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
        const std::vector<std::complex<double>>&) override;
};

#endif