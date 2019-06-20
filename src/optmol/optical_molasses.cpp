#include "optical_molasses.hpp"

const std::string DEFAULT_CFG_FILE = "config/params_optmol.cfg";
const std::string OUTPUT_DIR = "output/optmol";
const std::string ENERGY_OUTFILEBASE = "avgKE.out";
const std::string SPEED_DISTR_OUTFILEBASE = "speed_distr.out";
const unsigned OUTFILENAME_PRECISION = 3;

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cout << "Usage: ./optical_molasses <particle species> [<config file>]"
            << std::endl;
        return 1;
    }
    // Read in a possible config file
    std::string cfg_file = DEFAULT_CFG_FILE;
    if(argc > 2) {
        cfg_file = std::string(argv[2]);
    }

    // Get params
    // Accepted particle species strings:
    // "BePlus"
    // "Rb"
    PhysicalParams params(std::string(argv[1]), cfg_file);
    
    // Print out the params
    params.print();

    // Thermal velocity standard deviation of velocity components is sqrt(kT/m)
    double thermal_v_stddev = sqrt(K_BOLTZMANN*params.initial_temp/params.mass);

    // Initialize randomizer object with generator, thermal stddev,
    // and particle number
    // Set up RNG
    pcg32 generator(pcg_extras::seed_seq_from<std::random_device>{});
    // std::mt19937 generator(std::random_device{}());
    RandProcesses<pcg32> rng(generator, thermal_v_stddev, params.n_particles);
    
    // Initialize velocities to thermal distribution
    std::vector< std::vector<double> > v_particles(params.n_particles,
        std::vector<double>(3));
    for(auto i = v_particles.begin(); i != v_particles.end(); ++i) {
        for(auto vi = i->begin(); vi != i->end(); ++vi) {
            *vi = rng.rand_thermal_velocity();
        }
    }
    /// For output consistency with a single particle, force to have exactly
    /// the thermal energy
    if(params.n_particles == 1) {
        v_particles.back()[0] = thermal_v_stddev;
        v_particles.back()[1] = thermal_v_stddev;
        v_particles.back()[2] = thermal_v_stddev;
    }
    ///

    // Output files
    std::ostringstream suffix_ss;
    suffix_ss << std::setprecision(OUTFILENAME_PRECISION)
        << "N" << params.n_particles
        << "_Density" << params.particle_density
        << "_Omega" << params.rabi_freq_per_decay_rate
        << "_Delta" << params.initial_detuning_per_decay_rate << "to"
        << params.final_detuning_per_decay_rate
        << "_RampRate" << params.detuning_ramp_rate_natl_units
        << "_Temp" << params.initial_temp;

    std::ofstream energy_outfile(fullfile(tag_filename(
        ENERGY_OUTFILEBASE, suffix_ss.str(), params.particle_species),
        OUTPUT_DIR));
    // Number of energy snapshots to take
    unsigned n_snapshots = 1001;
    unsigned steps_between_snapshots = params.n_time_steps / (n_snapshots - 1);
    // Initial average kinetic energy
    energy_outfile << 0 << " "
        << calc_avg_kinetic_energy(v_particles, params.mass)/K_BOLTZMANN;

    // Initial speed distribution
    std::ofstream speed_init_outfile(fullfile(tag_filename(
        SPEED_DISTR_OUTFILEBASE, {"initial", suffix_ss.str()},
        params.particle_species),
        OUTPUT_DIR
    ));
    for(auto vp: v_particles) {
        speed_init_outfile << sqrt(sqr(vp[0]) + sqr(vp[1]) + sqr(vp[2]))
            << " ";
    }

    /// TIMING
    auto start = std::chrono::system_clock::now();
    ///

    /// COUNTING
    unsigned n_heat = 0, n_cool = 0;
    unsigned long n_collisions = 0;
    ///

    // Run over each time step
    for(unsigned i = 0; i < params.n_time_steps; ++i) {
        // Ramped detuning and photon wavenumber
        double detuning = calc_ramp((i+1)*params.dt,
            params.initial_detuning, params.final_detuning,
            params.detuning_ramp_rate);
        double laser_wavenumber = PhysicalParams::calc_laser_wavenumber(
            params.resonant_wavenumber, detuning);
        // velocity kick from a single photon absorption/emission
        double v_kick = HBAR*laser_wavenumber / params.mass;

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
                double absorb_rate = PhysicalParams::calc_absorb_rate(
                    params.decay_rate, params.rabi_freq, doppler_detuning);

                // Decide whether or not to absorb a photon
                if(!rng.rand_success_with_prob(
                    absorb_rate*params.dt - sqr(absorb_rate*params.dt)/2)) {
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
                double cos_theta, phi;
                std::tie(cos_theta, phi) = rng.rand_dir();
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
        // Note that scattering particles conserves kinetic energy,
        // so it doesn't have to be recomputed later
        double avgKE = calc_avg_kinetic_energy(v_particles, params.mass);
        if((i+1) % steps_between_snapshots == 0) {
            energy_outfile << std::endl << (i+1)*params.dt << " " <<
                avgKE/K_BOLTZMANN;
        }

        // Scatter some number of particles if possible
        if(params.n_particles > 1 && params.scatter_coeff > 0) {
            for(unsigned i_scat = 0; i_scat < params.collisions_per_step; ++i_scat) {
                // Choose two particles to scatter
                auto idxs = rng.rand_idx_pair();

                // Decide whether to scatter or not
                double rel_speed = calc_rel_speed(
                    v_particles[idxs.first], v_particles[idxs.second]);
                // 1+ to keep the argument above 1
                double coulomb_log = log(1 + 12*M_PI/cube(ELEMENTARY_CHARGE)
                    * sqrt(8*cube(VACUUM_PERMITTIVITY*avgKE)/(27*params.particle_density)));
                double scatter_prob = params.scatter_coeff*coulomb_log*params.dt
                    /cube(rel_speed);

                if(!rng.rand_success_with_prob(scatter_prob)) {
                    continue;
                }
                n_collisions += 1;

                // Carry on with scattering the pair
                // Get a random direction for scattering
                double cos_theta, phi;
                std::tie(cos_theta, phi) = rng.rand_dir();
                double sin_theta = sqrt(1 - sqr(cos_theta));
                
                auto scattered_vels = scatter_pair(
                    v_particles[idxs.first], v_particles[idxs.second],
                    {sin_theta*cos(phi), sin_theta*sin(phi), cos_theta}
                    );
                v_particles[idxs.first] = scattered_vels.first;
                v_particles[idxs.second] = scattered_vels.second;
            }
        }
    }
    energy_outfile.close();

    // Final speed distribution
    std::ofstream speed_final_outfile(fullfile(tag_filename(
        SPEED_DISTR_OUTFILEBASE, {"final", suffix_ss.str()},
        params.particle_species),
        OUTPUT_DIR
    ));
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
    std::cout << "Average collision success rate per time step: "
        << static_cast<double>(n_collisions)
            /(params.collisions_per_step*params.n_time_steps)
        << std::endl;
    std::cout << "Collision rate/max absorption rate: "
        << n_collisions / (params.duration_by_max_absorb_rate) << std::endl;
    ///
}

double calc_ramp(double t, double init, double final, double rate) {
    return std::max(std::min(init, final),
        std::min(std::max(init, final),
            init + rate*t
        ));
}

double calc_avg_kinetic_energy(const std::vector< std::vector<double> >& velocities,
    double mass) {
    double sumKE = 0;
    for(auto v: velocities) {
        sumKE += 0.5*mass*(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
    }
    return sumKE / velocities.size();
}

double calc_rel_speed(
    const std::vector<double>& v1, const std::vector<double>& v2) {
    return sqrt(sqr(v1[0]-v2[0]) + sqr(v1[1]-v2[1]) + sqr(v1[2]-v2[2]));
}

std::pair< std::vector<double>, std::vector<double> > scatter_pair(
    const std::vector<double>& v1, const std::vector<double>& v2,
    const std::vector<double>& rand_dir) {
    double rel_speed = calc_rel_speed(v1, v2);
    
    std::vector<double> new_v1, new_v2;
    new_v1.reserve(v1.size());
    new_v2.reserve(v2.size());
    for(unsigned i = 0; i < v1.size(); ++i) {
        new_v1.push_back((v1[i] + v2[i] + rel_speed*rand_dir[i]) / 2);
        new_v2.push_back((v1[i] + v2[i] - rel_speed*rand_dir[i]) / 2);
    }
    return std::make_pair(new_v1, new_v2);
}