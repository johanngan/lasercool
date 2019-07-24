#ifndef HSAWTOOTHPSI_HPP_
#define HSAWTOOTHPSI_HPP_

#include <string>
#include <vector>
#include <complex>
#include "readcfg.hpp"

// Two-level Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency
struct HSawtoothPsi {
    double rabi_freq_per_decay, detun_amp_per_decay, detun_freq_per_decay;
    double transition_angfreq_per_decay, omega1_per_decay, omega2_per_decay;

    HSawtoothPsi(std::string);
    double cumulative_phase(double);
    // Derivative operator to be passed to the timestepper
    std::vector<std::complex<double>> operator()(double,
        const std::vector<std::complex<double>>&);
};

#endif