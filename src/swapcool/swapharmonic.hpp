#ifndef SWAP_HARMONIC_HPP_
#define SWAP_HARMONIC_HPP_

#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <complex>
#include <vector>
#include <chrono>
#include "HHarmonic.hpp"
#include "DensMatHandler.hpp"
#include "lasercool/readcfg.hpp"
#include "lasercool/iotag.hpp"
#include "lasercool/timestepping.hpp"
#include "lasercool/fundconst.hpp"

// Generate a thermal state
std::vector<std::complex<double>> thermal_state(double, const HHarmonic&);
// Print out information about the system
void print_system_info(const HHarmonic&, double, double, bool, double, double);
// Calculate the average h value of a state
double calc_havg(const std::vector<std::complex<double>>&,
    const DensMatHandler&);
// Calculate the average h value within the population that hasn't leaked yet
double calc_havg_unleaked(const std::vector<std::complex<double>>&,
    const DensMatHandler&);
// Write state info to a file given the density matrix at a fixed time
void write_state_info(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const DensMatHandler&);
// Write the h-distribution at a fixed time to a file in tall format
void write_hdist(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const DensMatHandler&);

#endif