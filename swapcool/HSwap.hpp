#ifndef HSWAP_HPP_
#define HSWAP_HPP_

#include <string>
#include <vector>
#include <complex>
#include "read_config.hpp"

// Abstract class.
// Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency, and an exponential soft switch for the
// Rabi frequency.
struct HSwap{
    const double HBAR;
    double branching_ratio; // To the "low" (but not ground) state
    double enable_decay;    // 1 for enabled and 0 for disabled
    double rabi_freq_per_decay, detun_amp_per_decay, detun_freq_per_decay;
    double rabi_switch_coeff, rabi_switch_power;
    double transition_angfreq_per_decay;

    HSwap(std::string);

    // Rabi frequency soft switch on/off at a given (decay rate)*time
    // starting from 0 at gamma*t = 0 (mod gamma/f)
    double rabi_softswitch(double) const;
    // Detuning with sawtooth oscillation at a given (decay rate)*time,
    // starting from the minimum value at gamma*t = 0
    double detun_per_decay(double) const;
    // After some time (decay rate)*t
    // Assumes detuning chirp frequency is nonzero
    double cumulative_phase(double) const;

    // Transforms the coefficients solved for in the rotating wave
    // approximation back to the actual density matrix values;
    // i.e. put the oscillation back in.
    virtual std::vector<std::complex<double>> density_matrix(
        double, const std::vector<std::complex<double>>&) const = 0;

    // Derivative operator to be passed to the timestepper
    virtual std::vector<std::complex<double>> operator()(double,
        const std::vector<std::complex<double>>&) const = 0;
};

#endif