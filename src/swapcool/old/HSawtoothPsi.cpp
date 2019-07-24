#include "HSawtoothPsi.hpp"
using namespace std::complex_literals;
const double HBAR = 1.054571800e-34;

HSawtoothPsi::HSawtoothPsi(std::string fname) {
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
    omega1_per_decay = low_energy/(HBAR*decay_rate);
    omega2_per_decay = high_energy/(HBAR*decay_rate);
    transition_angfreq_per_decay = omega2_per_decay - omega1_per_decay;
}

// After some time (decay rate)*t
// Assumes detuning chirp frequency is nonzero
double HSawtoothPsi::cumulative_phase(double gt) {
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
std::vector<std::complex<double>> HSawtoothPsi::operator()(double gt,
    const std::vector<std::complex<double>>& psi) {
    double rabi_cos = rabi_freq_per_decay*cos(cumulative_phase(gt));
    return {
        -1i*(omega1_per_decay*psi[0] + rabi_cos*psi[1]),
        -1i*(rabi_cos*psi[0] + omega2_per_decay*psi[1])
    };
}
