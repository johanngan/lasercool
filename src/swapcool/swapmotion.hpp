#ifndef SWAPMOTION_HPP_
#define SWAPMOTION_HPP_

#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <complex>
#include <vector>
#include "HMotion.hpp"
#include "lasercool/readcfg.hpp"
#include "lasercool/iotag.hpp"
#include "lasercool/timestepping.hpp"

// Write state info to a file given the density matrix at a fixed time
void write_state_info(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const HMotion&);
// Write the k-distribution at a fixed time to a file in tall format
void write_kdist(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const HMotion&);

#endif