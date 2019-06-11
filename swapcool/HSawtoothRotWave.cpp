#include "HSawtoothRotWave.hpp"
using namespace std::complex_literals;
const double HBAR = 1.054571800e-34;

HSawtoothRotWave::HSawtoothRotWave(std::string fname) {
    double decay_rate, low_energy, high_energy;
    load_params(fname,
        {
            {"spontaneous_decay_rate", &decay_rate},
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

double HSawtoothRotWave::rabi_softswitch(double gt) {
    double _;
    return rabi_freq_per_decay*std::exp(-rabi_switch_coeff*std::pow(
        std::abs(2*modf(detun_freq_per_decay*gt, &_) - 1), rabi_switch_power));
}

double HSawtoothRotWave::detun_per_decay(double gt) {
    double _;
    return detun_amp_per_decay * (2*modf(detun_freq_per_decay*gt, &_) - 1);
}

double HSawtoothRotWave::cumulative_phase(double gt) {
    double ncycles;
    double cycle_completion = modf(gt*detun_freq_per_decay, &ncycles);
    // Phase from full cycles + the phase from the current one
    return (ncycles * transition_angfreq_per_decay  // Full cycles
        + cycle_completion  // Current cycle
            * (transition_angfreq_per_decay
                + detun_amp_per_decay*(cycle_completion - 1))
        ) / detun_freq_per_decay;
}

std::vector<std::complex<double>> HSawtoothRotWave::density_matrix(
    double gt, const std::vector<std::complex<double>>& coefficients) {
    std::complex<double> cexp = std::exp(1i*cumulative_phase(gt));
    return {
        coefficients[0],
        coefficients[1]*cexp,
        coefficients[2]*std::conj(cexp),
        coefficients[3]
    };
}

// Assumes row major format
std::vector<std::complex<double>> HSawtoothRotWave::operator()(double gt,
    const std::vector<std::complex<double>>& rho_c) {
    // 1/(i*HBAR) * [H, rho_c] + L(rho_c) from the master equation
    return {
        rho_c[3] + 0.5i*rabi_softswitch(gt)*(rho_c[1]-rho_c[2]),
        -0.5*rho_c[1] + 1i*(0.5*rabi_softswitch(gt)*(rho_c[0]-rho_c[3])
            - detun_per_decay(gt)*rho_c[1]),
        -0.5*rho_c[2] - 1i*(0.5*rabi_softswitch(gt)*(rho_c[0]-rho_c[3])
            - detun_per_decay(gt)*rho_c[2]),
        -rho_c[3] - 0.5i*rabi_softswitch(gt)*(rho_c[1]-rho_c[2])
    };
}
