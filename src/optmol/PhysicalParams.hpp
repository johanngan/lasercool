#ifndef PHYSICALPARAMS_HPP_
#define PHYSICALPARAMS_HPP_

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include "lasercool/readcfg.hpp"
#include "constants.hpp"
#include "mathutil.hpp"

// Read, calculate, and hold relevant physical parameters
struct PhysicalParams {
    // Particle species stuff
    std::string particle_species;
    double mass, decay_rate, resonant_wavenumber;

    // Config file stuff
    double rabi_freq_per_decay_rate, initial_detuning_per_decay_rate, 
        final_detuning_per_decay_rate, detuning_ramp_rate_natl_units;    
    double initial_temp;
    double dt_by_max_absorb_rate, duration_by_max_absorb_rate;
    unsigned n_particles;
    double particle_density;

    // Stuff in SI units
    double rabi_freq, initial_detuning, final_detuning, detuning_ramp_rate;
    double dt, duration;

    // Calculated stuff
    double max_absorb_rate;
    unsigned n_time_steps;
    unsigned collisions_per_step;
    double scatter_coeff;

    // Initialize with a given particle species and
    // read other parameters in from a config file
    PhysicalParams(std::string, std::string);
    // Print out params to console
    void print();

    // Photon absorption rate for given detuning
    inline static double calc_absorb_rate(double decay_rate, double rabi_freq,
        double detuning=0) {
        return 0.25*decay_rate*sqr(rabi_freq)
            / (sqr(detuning) + 0.5*sqr(rabi_freq) + 0.25*sqr(decay_rate));
    }

    inline static double calc_laser_wavenumber(double resonant_wavenumber,
        double detuning) {
        return resonant_wavenumber + detuning/SPEED_OF_LIGHT;
    }

    // Equilibrium temperature, nonzero due to shot heating
    inline static double expected_min_temp(double decay_rate,
        double detuning) {
        return -0.125*HBAR/K_BOLTZMANN
            * (sqr(decay_rate) + 4*sqr(detuning))/detuning;
    }

    // Optimal detuning in the high temperature limit, which leads to the
    // greatest instantaneous cooling rate at a given temperature
    inline static double optimal_detuning(double temp, double mass,
        double wavenumber) {
        return -wavenumber * sqrt(K_BOLTZMANN*temp/mass);
    }
};

#endif