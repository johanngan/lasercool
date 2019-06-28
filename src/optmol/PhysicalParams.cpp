#include "PhysicalParams.hpp"

PhysicalParams::PhysicalParams(std::string species, std::string fname):
    particle_species(species) {
    // Set particle-species–specific parameters
    // Resonant wavenumber is the wavenumber of the atomic transition
    if(particle_species == "BePlus") {
        mass = MASS_BE_PLUS;
        decay_rate = DECAY_RATE_BE_PLUS;
        resonant_wavenumber = WAVENUMBER_BE_2S2P;
    } else if(particle_species == "Rb") {
        mass = MASS_RB;
        decay_rate = DECAY_RATE_RB;
        resonant_wavenumber = WAVENUMBER_RB;
    } else {
        throw std::invalid_argument("Invalid particle species");
    }

    // Read in parameters from cfg file
    double n_particles_double;
    load_params(fname,
        {
            {"rabi_frequency", &rabi_freq_per_decay_rate},
            {"initial_detuning", &initial_detuning_per_decay_rate},
            {"final_detuning", &final_detuning_per_decay_rate},
            {"detuning_ramp_rate", &detuning_ramp_rate_natl_units},
            {"initial_temperature", &initial_temp},
            {"time_step", &dt_by_max_absorb_rate},
            {"duration", &duration_by_max_absorb_rate},
            {"n_particles", &n_particles_double},
            {"particle_density", &particle_density}
        }
    );
    n_particles = static_cast<unsigned>(n_particles_double);
    rabi_freq = rabi_freq_per_decay_rate * decay_rate;
    // Set time scale in terms of maximum photon absorption rate
    max_absorb_rate = calc_absorb_rate(decay_rate, rabi_freq);
    dt = dt_by_max_absorb_rate / max_absorb_rate;
    duration = duration_by_max_absorb_rate / max_absorb_rate;
    n_time_steps = static_cast<unsigned>(ceil(
        duration_by_max_absorb_rate / dt_by_max_absorb_rate));
    
    // Precompute certain values for the scattering rate
    // So each particle has on average one collision per time step
    collisions_per_step = n_particles / 2;
    // Unchanging coefficient on the scattering probability, where
    // P(scatter in interval dt) = coeff*(coulomb log)*dt/(relative speed)^3
    // See http://www.physics.purdue.edu/~robichf/papers/PoP10_2217.pdf
    scatter_coeff = particle_density*pow(ELEMENTARY_CHARGE, 4)
        / (M_PI*sqr(VACUUM_PERMITTIVITY*mass));

    // Defaults and conversion to SI //
    if(std::isnan(final_detuning_per_decay_rate)) {
        final_detuning_per_decay_rate = -0.5;   // Gives Doppler temperature
    }
    final_detuning = final_detuning_per_decay_rate * decay_rate;
    if(std::isnan(initial_detuning_per_decay_rate)) {
        // Gives highest cooling rate at high temperatures
        initial_detuning_per_decay_rate = optimal_detuning(initial_temp, mass,
            calc_laser_wavenumber(resonant_wavenumber, final_detuning))
            / decay_rate;
    }
    initial_detuning = initial_detuning_per_decay_rate * decay_rate;
    if(std::isnan(detuning_ramp_rate_natl_units)) {
        // Ramp the entire time
        detuning_ramp_rate_natl_units =
            (final_detuning_per_decay_rate - initial_detuning_per_decay_rate)
            / duration_by_max_absorb_rate;
    }
    detuning_ramp_rate = detuning_ramp_rate_natl_units
        * decay_rate * max_absorb_rate;
}

void PhysicalParams::print() {
    // Output parameters
    std::cout
        << "Parameters:" << std::endl
        << "    Particle species: " << particle_species << std::endl
        << "    N: " << n_particles << std::endl
        << "    Density: " << particle_density << std::endl
        << "    Rabi frequency per decay rate: " << rabi_freq_per_decay_rate
        << std::endl
        << "    Initial detuning per decay rate: "
        << initial_detuning_per_decay_rate << std::endl
        << "    Final detuning per decay rate: "
        << final_detuning_per_decay_rate << std::endl
        << "    Detuning ramp rate * max absorption rate / decay rate: "
        << detuning_ramp_rate_natl_units << std::endl
        << "    Final laser wavelength: "
        << 2*M_PI/calc_laser_wavenumber(resonant_wavenumber, final_detuning)
             * 1e9 << " nm" << std::endl
        << "    Temperature: " << initial_temp << " K" << std::endl
        << "    Time step * max absorption rate: " << dt_by_max_absorb_rate
        << std::endl
        << "    Duration * max absorption rate: " << duration_by_max_absorb_rate
        << std::endl;
    // Output useful, theoretically calculated quantities related to optimization
    std::cout << "Optimal initial detuning per decay rate: "
        << optimal_detuning(initial_temp, mass,
            calc_laser_wavenumber(resonant_wavenumber, final_detuning))
            / decay_rate
        << std::endl
        << "Expected optical molasses equilibrium temperature: "
        << expected_min_temp(decay_rate, final_detuning) << " K"
        << std::endl;
}
