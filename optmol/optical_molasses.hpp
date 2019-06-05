#ifndef OPTICAL_MOLASSES_HPP_
#define OPTICAL_MOLASSES_HPP_

#include <cmath>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include "constants.hpp"
#include "mathutil.hpp"
#include "PhysicalParams.hpp"
#include "RandProcesses.hpp"
#include "pcg_random.hpp"

// Calculate a ramped quantity over time given the initial and final values,
// and the ramp rate
double calc_ramp(double, double, double, double);
// Calculate average kinetic energy of an ensemble
double calc_avg_kinetic_energy(const std::vector< std::vector<double> >&,
    double);
// Compute the relative speed between two particles
double calc_rel_speed(const std::vector<double>&, const std::vector<double>&);
// Scatter two particles in a collision
std::pair< std::vector<double>, std::vector<double> > scatter_pair(
    const std::vector<double>&, const std::vector<double>&,
    const std::vector<double>&);

// Tag a file name with some suffix
std::string tag_filename(std::string, std::string,
    std::string prefix="", std::string separator="_");

#endif