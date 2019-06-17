#include "HSwap.hpp"

HSwap::HSwap(std::string fname):HBAR(1.054571800e-34) {
    double decay_rate, low_energy, high_energy;
    load_params(fname,
        {
            {"spontaneous_decay_rate", &decay_rate},
            {"enable_decay", &enable_decay},
            {"branching_ratio", &branching_ratio},
            {"low_energy_level", &low_energy},
            {"high_energy_level", &high_energy},
            {"rabi_frequency", &rabi_freq_per_decay},
            {"rabi_switch_coeff", &rabi_switch_coeff},
            {"rabi_switch_power", &rabi_switch_power},
            {"detuning_amplitude", &detun_amp_per_decay},
            {"detuning_frequency", &detun_freq_per_decay}
        }
    );
    transition_angfreq_per_decay = (high_energy - low_energy)/(HBAR*decay_rate);
}

double HSwap::rabi_softswitch(double gt) const {
    double _;
    return rabi_freq_per_decay*std::exp(-rabi_switch_coeff*std::pow(
        std::abs(2*modf(detun_freq_per_decay*gt, &_) - 1), rabi_switch_power));
}

double HSwap::detun_per_decay(double gt) const {
    double _;
    return detun_amp_per_decay * (2*modf(detun_freq_per_decay*gt, &_) - 1);
}

double HSwap::cumulative_phase(double gt) const {
    double ncycles;
    double cycle_completion = modf(gt*detun_freq_per_decay, &ncycles);
    // Phase from full cycles + the phase from the current one
    return (ncycles * transition_angfreq_per_decay  // Full cycles
        + cycle_completion  // Current cycle
            * (transition_angfreq_per_decay
                + detun_amp_per_decay*(cycle_completion - 1))
        ) / detun_freq_per_decay;
}
