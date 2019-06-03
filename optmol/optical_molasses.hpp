#ifndef OPTICAL_MOLASSES_HPP_
#define OPTICAL_MOLASSES_HPP_

#include <cmath>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include "constants.hpp"
#include "read_config.hpp"
#include "pcg_random.hpp"

// Square a number
double sqr(double);
// Calculate the absorption rate of photons
double calc_absorb_rate(double, double, double detuning=0);
// Calculate the laser light wavenumber
double calc_laser_wavenumber(double, double);
// Calculate average kinetic energy of an ensemble
double calc_avg_kinetic_energy(const std::vector< std::vector<double> >&,
    double);
// Calculate optical molasses theoretical equilibrium temperature
double expected_min_temp(double, double);
// Tag a file name with some suffix
std::string tag_filename(std::string, std::string,
    std::string prefix="", std::string separator="_");

#endif