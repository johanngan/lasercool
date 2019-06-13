#ifndef HSAWTOOTHROTWAVE_HPP_
#define HSAWTOOTHROTWAVE_HPP_

#include <string>
#include <vector>
#include <complex>
#include "read_config.hpp"

// Two-level Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency, under the rotating wave approximation
struct HSawtoothRotWave {
    const unsigned nstates; // "matrix dimension"
    double rabi_freq_per_decay, detun_amp_per_decay, detun_freq_per_decay;
    double rabi_switch_coeff, rabi_switch_power;
    double transition_angfreq_per_decay;

    HSawtoothRotWave(std::string);
    // Convert matrix subscripts to linear indexes (row-major format)
    unsigned subidx(unsigned, unsigned);

    // Rabi frequency soft switch on/off at a given (decay rate)*time
    // starting from 0 at gamma*t = 0 (mod gamma/f)
    double rabi_softswitch(double);
    // Detuning with sawtooth oscillation at a given (decay rate)*time,
    // starting from the minimum value at gamma*t = 0
    double detun_per_decay(double);
    // After some time (decay rate)*t
    // Assumes detuning chirp frequency is nonzero
    double cumulative_phase(double);
    // Transforms the coefficients solved for in the rotating wave
    // approximation back to the actual density matrix values;
    // i.e. put the oscillation back in.
    std::vector<std::complex<double>> density_matrix(
        double, const std::vector<std::complex<double>>&);

    // Derivative operator to be passed to the timestepper
    std::vector<std::complex<double>> operator()(double,
        const std::vector<std::complex<double>>&);
};

#endif