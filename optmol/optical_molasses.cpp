#include "optical_molasses.hpp"

const std::string CFG_FILE = "params.cfg";
const std::string OUTPUT_DIR = "output";
const std::string ENERGY_OUTFILEBASE = "avgKE.out";
const std::string SPEED_DISTR_OUTFILEBASE = "speed_distr.out";

int main(int argc, char** argv) {
    if(argc != 2) {
        std::cout << "Usage: ./optical_molasses <particle species>"
            << std::endl;
        return 1;
    }

    // Particle species string
    // "BePlus"
    // "Rb"
    std::string particle_species(argv[1]);

    // Resonant wavenumber is the wavenumber of the atomic transition
    double mass, decay_rate, resonant_wavenumber;
    if(particle_species == "BePlus") {
        mass = MASS_BE_PLUS;
        decay_rate = DECAY_RATE_BE_PLUS;
        resonant_wavenumber = WAVENUMBER_BE_2S2P;
    } else if(particle_species == "Rb") {
        mass = MASS_RB;
        decay_rate = DECAY_RATE_RB;
        resonant_wavenumber = WAVENUMBER_RB;
    } else {
        std::cout << "ERROR: Invalid particle species." << std::endl;
        return 1;
    }

    // Read in parameters
    double rabi_freq_per_decay_rate, detuning_per_decay_rate, initial_temp,
        dt, duration, n_particles_double;
    unsigned n_particles;
    load_params(CFG_FILE,
        {
            {"rabi_frequency", &rabi_freq_per_decay_rate},
            {"detuning", &detuning_per_decay_rate},
            {"initial_temperature", &initial_temp},
            {"time_step", &dt},
            {"duration", &duration},
            {"n_particles", &n_particles_double}
        }
    );
    n_particles = static_cast<unsigned>(n_particles_double);
    // Set Rabi frequency and detuning in terms of spontaneous decay rate
    double rabi_freq = rabi_freq_per_decay_rate * decay_rate;
    double detuning = detuning_per_decay_rate * decay_rate;

    // Calculate the wavenumber of the laser light
    double laser_wavenumber = calc_laser_wavenumber(resonant_wavenumber,
        detuning);
    // velocity kick from a single photon absorption/emission
    double v_kick = HBAR*laser_wavenumber / mass;

    std::cout
        << "Parameters:" << std::endl
        << "    Particle species: " << particle_species << std::endl
        << "    N: " << n_particles << std::endl
        << "    Rabi frequency per decay rate: " << rabi_freq_per_decay_rate
        << std::endl
        << "    Detuning per decay rate: " << detuning_per_decay_rate
        << std::endl
        << "    Laser wavelength: " << 2*M_PI/laser_wavenumber * 1e9 << " nm"
        << std::endl
        << "    Temperature: " << initial_temp << " K" << std::endl
        << "    Time step * max decay rate: " << dt << std::endl
        << "    Duration * max decay rate: " << duration << std::endl;
    std::cout << "Expected optical molasses equilibrium temperature: "
        << expected_min_temp(decay_rate, detuning) << " K"
        << std::endl;

    // Calculate number of time steps for speed
    unsigned n_time_steps = static_cast<unsigned>(ceil(duration / dt));

    // Set time step and duration in terms of the max absorption rate
    double max_absorb_rate = calc_absorb_rate(decay_rate, rabi_freq);
    dt /= max_absorb_rate;
    duration /= max_absorb_rate;

    // Set up RNG
    std::mt19937 generator(std::random_device{}());

    // Initialize particles to a thermal distribution
    // Gaussian in each velocity component, with mean 0 and
    // standard deviation sqrt(kT/m)
    double mean = 0, stddev = sqrt(K_BOLTZMANN*initial_temp/mass);
    std::normal_distribution<> normal_dist(mean, stddev);
    std::vector< std::vector<double> > v_particles(n_particles,
        std::vector<double>(3));
    for(auto i = v_particles.begin(); i != v_particles.end(); ++i) {
        for(auto vi = i->begin(); vi != i->end(); ++vi) {
            *vi = normal_dist(generator);
        }
    }

    /// For output consistency with a single particle, force to have exactly
    /// the thermal energy
    if(n_particles == 1) {
        v_particles.back()[0] = sqrt(K_BOLTZMANN*initial_temp/mass);
        v_particles.back()[1] = sqrt(K_BOLTZMANN*initial_temp/mass);
        v_particles.back()[2] = sqrt(K_BOLTZMANN*initial_temp/mass);
    }
    ///

    // Output files
    std::ostringstream suffix_ss;
    suffix_ss << "N" << n_particles
        << "_Omega" << rabi_freq_per_decay_rate
        << "_Delta" << detuning_per_decay_rate
        << "_Temp" << initial_temp;

    std::ofstream energy_outfile(OUTPUT_DIR + "/" + tag_filename(
        ENERGY_OUTFILEBASE, suffix_ss.str(), particle_species));
    // Number of energy snapshots to take
    unsigned n_snapshots = 1001;
    unsigned steps_between_snapshots = n_time_steps / (n_snapshots - 1);
    // Initial average kinetic energy
    energy_outfile << 0 << " "
        << calc_avg_kinetic_energy(v_particles, mass)/K_BOLTZMANN;

    // Initial speed distribution
    std::ofstream speed_init_outfile(OUTPUT_DIR + "/" + tag_filename(
        tag_filename(SPEED_DISTR_OUTFILEBASE, "initial"),
        suffix_ss.str(), particle_species));
    for(auto vp: v_particles) {
        speed_init_outfile << sqrt(sqr(vp[0]) + sqr(vp[1]) + sqr(vp[2]))
            << " ";
    }

    /// TIMING
    auto start = std::chrono::system_clock::now();
    ///

    /// COUNTING
    unsigned n_heat = 0, n_cool = 0;
    ///

    // [0, 1) uniform distribution
    std::uniform_real_distribution<> uniform_dist;
    // Run over each time step
    for(unsigned i = 0; i < n_time_steps; ++i) {
        // Run over each particle
        for(auto vp = v_particles.begin(); vp != v_particles.end(); ++vp) {
            // Iterate over each of the 6 lasers
            // Goes through -x, +x, -y, +y, -z, +z
            int direction = 1;
            for(unsigned j = 0; j < 6; ++j) {
                int component = j / 2;
                direction *= -1;

                // Doppler-shifted detuning
                double doppler_detuning = detuning
                    - direction*laser_wavenumber*(*vp)[component];
                double absorb_rate = calc_absorb_rate(decay_rate,
                    rabi_freq, doppler_detuning);

                // Decide whether or not to absorb a photon
                if(uniform_dist(generator) >= absorb_rate*dt) {
                    continue;
                }

                if((*vp)[component]*direction > 0) {
                    n_heat++;
                } else {
                    n_cool++;
                }

                // Photon absorbed //
                // Absorption kick
                (*vp)[component] += direction*v_kick;
                // Get a random direction for emission
                double cos_theta = 2*uniform_dist(generator) - 1;
                double phi = 2*M_PI*uniform_dist(generator);
                double sin_theta = sqrt(1 - sqr(cos_theta));
                
                // Emission kick
                (*vp)[0] += v_kick*sin_theta*cos(phi);
                (*vp)[1] += v_kick*sin_theta*sin(phi);
                (*vp)[2] += v_kick*cos_theta;
            }
            // Insert the desired measurement calculations //
        }
        // Insert the desired measurement calculations //

        // Average kinetic energy
        if((i+1) % steps_between_snapshots == 0) {
            energy_outfile << std::endl << (i+1)*dt << " " <<
                calc_avg_kinetic_energy(v_particles, mass)/K_BOLTZMANN;
        }
    }
    energy_outfile.close();

    // Final speed distribution
    std::ofstream speed_final_outfile(OUTPUT_DIR + "/" + tag_filename(
        tag_filename(SPEED_DISTR_OUTFILEBASE, "final"),
        suffix_ss.str(), particle_species));
    for(auto vp: v_particles) {
        speed_final_outfile << sqrt(sqr(vp[0]) + sqr(vp[1]) + sqr(vp[2]))
            << " ";
    }

    ///
    std::chrono::duration<double> total_seconds =
        std::chrono::system_clock::now() - start;
    std::cout << "Total runtime: " << total_seconds.count() << " s" << std::endl;
    std::cout << "Number of heating events: " << n_heat << std::endl
        << "Number of cooling events: " << n_cool << std::endl;
    ///
}

double sqr(double x) {
    return x*x;
}

double calc_absorb_rate(double decay_rate, double rabi_freq, double detuning) {
    return 0.25*decay_rate*sqr(rabi_freq)
        / (sqr(detuning) + 0.5*sqr(rabi_freq) + 0.25*sqr(decay_rate));
}

double calc_laser_wavenumber(double resonant_wavenumber, double detuning) {
    return resonant_wavenumber + detuning/SPEED_OF_LIGHT;
}

double calc_avg_kinetic_energy(const std::vector< std::vector<double> >& velocities,
    double mass) {
    double sumKE = 0;
    for(auto v: velocities) {
        sumKE += 0.5*mass*(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
    }
    return sumKE / velocities.size();
}

double expected_min_temp(double decay_rate, double detuning) {
    return -0.125*HBAR/K_BOLTZMANN * (sqr(decay_rate) + 4*sqr(detuning))/detuning;
}

std::string tag_filename(std::string filename, std::string suffix,
    std::string prefix, std::string separator) {
    std::string tagged = filename.insert(filename.rfind("."),
        separator + suffix);
    if(!prefix.empty()) {
        tagged = prefix + separator + tagged;
    }
    return tagged;
}