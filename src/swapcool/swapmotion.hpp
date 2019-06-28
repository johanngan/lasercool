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
#include "DensMatHandler.hpp"
#include "lasercool/readcfg.hpp"
#include "lasercool/iotag.hpp"
#include "lasercool/timestepping.hpp"

// Quality metric evaluation string
std::string evaluate_quality_metric(
    double,
    double low_thresh=1,
    double very_low_thresh=0.5,
    std::string okay_str="",
    std::string low_str="*",
    std::string very_low_str="**");
// Calculate the RMS k value of a state
double calc_krms(const std::vector<std::complex<double>>&,
    const DensMatHandler&);
// Modify the density matrix in preparation for a new cycle
void initialize_cycle(std::vector<std::complex<double>>&, const HMotion&);
// Write state info to a file given the density matrix at a fixed time
void write_state_info(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const DensMatHandler&);
// Write the k-distribution at a fixed time to a file in tall format
void write_kdist(std::ofstream&, double,
    const std::vector<std::complex<double>>&, const DensMatHandler&);

#endif