#include "HSawtooth.hpp"
using namespace std::complex_literals;
const double HBAR = 1.054571800e-34;

HSawtooth::HSawtooth(std::string fname) {
    double decay_rate, low_energy, high_energy;
    load_params(fname,
        {
            {"spontaneous_decay_rate", &decay_rate},
            {"low_energy_level", &low_energy},
            {"high_energy_level", &high_energy},
            {"rabi_frequency", &rabi_freq_per_decay},
            {"detuning_amplitude", &detun_amp_per_decay},
            {"detuning_frequency", &detun_freq_per_decay}
        }
    );
    transition_angfreq_per_decay = (high_energy - low_energy)/(HBAR*decay_rate);
}

// After some time (decay rate)*t
// Assumes detuning chirp frequency is nonzero
double HSawtooth::cumulative_phase(double gt) {
    double ncycles;
    double cycle_completion = modf(gt*detun_freq_per_decay, &ncycles);
    // Phase from full cycles + the phase from the current one
    return (ncycles * transition_angfreq_per_decay  // Full cycles
        + cycle_completion  // Current cycle
            * (transition_angfreq_per_decay
                + detun_amp_per_decay*(cycle_completion - 1))
        ) / detun_freq_per_decay;
}

// Assumes row major format
std::vector<std::complex<double>> HSawtooth::operator()(double gt,
    const std::vector<std::complex<double>>& rho) {
    double rabi_cos = rabi_freq_per_decay*cos(cumulative_phase(gt));
    // 1/(i*HBAR) * [H, rho] + L(rho) from the master equation
    return {
        rho[3] + 1i*(rho[1]-rho[2])*rabi_cos,
        -0.5*rho[1] + 1i*(
            transition_angfreq_per_decay*rho[1] + (rho[0]-rho[3])*rabi_cos),
        -0.5*rho[2] - 1i*(
            transition_angfreq_per_decay*rho[2] + (rho[0]-rho[3])*rabi_cos),
        -rho[3] - 1i*(rho[1]-rho[2])*rabi_cos
    };
}
