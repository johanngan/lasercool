#ifndef DENSMATHANDLER_HPP_
#define DENSMATHANDLER_HPP_

#include <string>
#include <complex>
#include <vector>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <exception>
#include "lasercool/readcfg.hpp"

// Handler for dealing with an efficiently stored density matrix for the
// 3-state SWAP system that takes advantage of hermiticity and incoherence
// between the "sink" state and the other states.
struct DensMatHandler {
    
    
    unsigned nint;  // number of internal states
    int kmax, kmin;   // range of tracked k values
    unsigned kstates;   // number of k states
    // Indexes of matrix elements that are actually stored
    // The key is the linear index as determined by subidx()
    // The value is the index of the actual density matrix vector instance
    std::unordered_map<unsigned, unsigned> idxmap;
    // Contains the list of matrix elements at subscript (nl, kl, nr, kr)
    // that are actually stored
    std::vector<std::tuple<unsigned, int, unsigned, int>> idxlist;

    DensMatHandler(std::string);

    // Convert state subscripts to linear indexes in the density matrix,
    // enumerated as |n-left, k-left><n-right, k-right|
    inline unsigned subidx(unsigned, int, unsigned, int) const;

    // Checks if an element at some subscript is stored
    bool has(unsigned, int, unsigned, int) const;

    // Get the matrix element at some subscript
    std::complex<double> ele(const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int) const;
    // Get a reference to the element at some subscript
    std::complex<double>& at(std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int) const;
    const std::complex<double>& at(const std::vector<std::complex<double>>&,
        unsigned, int, unsigned, int) const;

    // Total trace
    std::complex<double> totaltr(
        const std::vector<std::complex<double>>&) const;

    // Partial trace over k for a fixed n
    std::complex<double> partialtr_k(
        const std::vector<std::complex<double>>&, unsigned) const;

    // Partial trace over n for a fixed k
    std::complex<double> partialtr_n(
        const std::vector<std::complex<double>>&, int) const;
    
    // Trace of rho^2
    std::complex<double> purity(const std::vector<std::complex<double>>&) const;
};

#endif