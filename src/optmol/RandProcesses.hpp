// Collection of random physical processes, reusing a single RNG and
// distribution instances to save on the expensive initialization
#ifndef RANDPROCESSES_HPP_
#define RANDPROCESSES_HPP_

#include <cmath>
#include <random>
#include <utility>

template<typename rngtype>
class RandProcesses {
    private:
        rngtype generator;
        // Mean 0 and stddev given at initialization
        std::normal_distribution<> thermal_v_dist;
        // Uniform [0, 1)
        std::uniform_real_distribution<> uniform_dist;
        // Uniform [0, 2*pi). Saves extra multiplication from the U(0, 1)
        std::uniform_real_distribution<> uniform_phi_dist;
        // Uniform [-1, 1).
        std::uniform_real_distribution<> uniform_costheta_dist;
        // Uniform [0, n-1], n given at initialization
        std::uniform_int_distribution<> idx_dist;
    public:
        // Initialize with relevant physical parameters
        RandProcesses(rngtype generator, double thermal_v_stddev, unsigned n_particles):
            generator(generator), thermal_v_dist(0., thermal_v_stddev),
            uniform_dist(), uniform_phi_dist(0., 2*M_PI),
            uniform_costheta_dist(-1., 1.),
            idx_dist(0, n_particles-1) {}

        // Generate velocities from a thermal distribution
        double rand_thermal_velocity() {return thermal_v_dist(generator);}

        // Generate pair of two distinct indexes. Assumes n_particles > 1,
        // otherwise there'll be trouble
        std::pair<unsigned, unsigned> rand_idx_pair() {
            unsigned idx1 = idx_dist(generator), idx2;
            do {
                idx2 = idx_dist(generator);
            } while(idx1 == idx2);
            return std::make_pair(idx1, idx2);
        }

        // Rolls a random success or failure with given success probability
        bool rand_success_with_prob(double prob) {
            return uniform_dist(generator) < prob;
        }

        // Random direction on the unit sphere
        // returns a random (cos(theta) ~ U(-1, 1), phi ~ U(0, 2*pi))
        std::pair<double, double> rand_dir() {
            return std::make_pair(
                uniform_costheta_dist(generator), uniform_phi_dist(generator));
        }
};

#endif