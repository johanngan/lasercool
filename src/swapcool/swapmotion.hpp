#ifndef SWAPMOTION_HPP_
#define SWAPMOTION_HPP_

#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <complex>
#include <vector>
#include <chrono>
#include "HMotion.hpp"
#include "lasercool/readcfg.hpp"
#include "lasercool/iotag.hpp"
#include "lasercool/timestepping.hpp"

// Modify the density matrix in preparation for a new cycle
void initialize_cycle(std::vector<std::complex<double>>&, const HMotion&);
// Write state info to a file given the density matrix at a fixed time
void write_state_info(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const HMotion&);
// Write the k-distribution at a fixed time to a file in tall format
void write_kdist(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const HMotion&);

#endif