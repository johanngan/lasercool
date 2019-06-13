#ifndef HINT_HPP_
#define HINT_HPP_

#include "HSwap.hpp"

// Hamiltonian for sawtooth laser frequency oscillating about
// some transition frequency, under the rotating wave approximation,
// only paying attention to internal states
struct HInt : public HSwap {
    const unsigned nstates; // "matrix dimension"

    HInt(std::string fname):HSwap(fname), nstates(3) {}
    // Convert matrix subscripts to linear indexes (row-major format)
    unsigned subidx(unsigned, unsigned);

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