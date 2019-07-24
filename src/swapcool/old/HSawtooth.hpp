#ifndef HSAWTOOTH_HPP_
#define HSAWTOOTH_HPP_

#include <string>
#include <vector>
#include <complex>
#include "readcfg.hpp"

// Two-level Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency
struct HSawtooth {
    double rabi_freq_per_decay, detun_amp_per_decay, detun_freq_per_decay;
    double transition_angfreq_per_decay;

    HSawtooth(std::string);
    double cumulative_phase(double);
    // Derivative operator to be passed to the timestepper
    std::vector<std::complex<double>> operator()(double,
        const std::vector<std::complex<double>>&);
};

#endif