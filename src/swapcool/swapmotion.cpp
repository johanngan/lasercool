// SWAP laser setup considering both internal states and momentum states.
// The density matrix is stored in row major format, enumerated
// as |nl, kl><nr, kr|
#include "swapmotion.hpp"

const std::string DEFAULT_CFG_FILE = "config/params_swapcool.cfg";
const std::string OUTPUT_DIR = "output/swapcool/swapmotion";
const std::string RHO_OUTFILEBASE = "rho.out";
const std::string CYCLE_OUTFILEBASE = "cycles.out";
const std::string KDIST_OUTFILEBASE = "kdist.out";
const std::string KDIST_FINAL_OUTFILEBASE = "kdist_final.out";
const unsigned OUTFILENAME_PRECISION = 3;

int main(int argc, char** argv) {
    // Read in a possible config file
    std::string cfg_file = DEFAULT_CFG_FILE;
    if(argc > 1) {
        cfg_file = std::string(argv[1]);
    }

    double decay_rate, duration_by_decay, tol, init_k_double;
    load_params(cfg_file,
        {
            {"spontaneous_decay_rate", &decay_rate},
            {"duration", &duration_by_decay},
            {"tolerance", &tol},
            {"initial_momentum", &init_k_double}
        }
    );
    int init_k = static_cast<int>(init_k_double);

    // Form the derivative operator, in natural units
    // d(rho)/d(Gamma*t)
    HMotion hamil(cfg_file);

    // Form output files
    std::ostringstream oftag_ss;
    oftag_ss << std::setprecision(OUTFILENAME_PRECISION)
        << "A" << hamil.detun_amp_per_decay
        << "_f" << hamil.detun_freq_per_decay
        << "_Omega" << hamil.rabi_freq_per_decay
        << "_k" << init_k;
    std::ofstream cyclesout(fullfile(tag_filename(
        CYCLE_OUTFILEBASE, oftag_ss.str()),
        OUTPUT_DIR
    ));
    for(double gt = 0; gt < duration_by_decay;
        gt += 1./(100*hamil.detun_freq_per_decay)) {
        cyclesout << gt << " " << hamil.rabi_softswitch(gt) << " "
        << hamil.detun_per_decay(gt) << " " << hamil.cumulative_phase(gt)
        << std::endl;
    }
    cyclesout.close();

    std::ofstream rho_out(fullfile(tag_filename(
        RHO_OUTFILEBASE, oftag_ss.str()),
        OUTPUT_DIR
    ));
    std::ofstream kdistout(fullfile(tag_filename(
        KDIST_OUTFILEBASE, oftag_ss.str()),
        OUTPUT_DIR
    ));
    std::ofstream kdistfinalout(fullfile(tag_filename(
        KDIST_FINAL_OUTFILEBASE, oftag_ss.str()),
        OUTPUT_DIR
    ));

    // Solve the system
    std::vector<std::complex<double>> rho_c0(hamil.nmat());
    // Initialize all in one k-state
    rho_c0[hamil.subidx(1, init_k, 1, init_k)] = 1;

    // Solve the system in natural units with an adaptive RK method
    auto rho_c_solution = timestepping::odesolve(hamil, rho_c0,
        duration_by_decay, timestepping::AdaptiveRK(tol));

    // Write the solution in SI units
    // Write table header
    rho_out << "t |rho11| |rho22| |rho33| tr(rho) tr(rho^2)";
    for(int k = 0; k <= hamil.kmax; ++k) {
        rho_out << " |k" << k << "|";
    }
    rho_out << " |k_rms|";
    rho_out << std::endl;
    // Write table header
    std::string kdist_header = "t k P(k) P(n = 0, k), P(n = 1, k), P(n = 2, k)";
    kdistout << kdist_header << std::endl;
    for(auto point: rho_c_solution) {
        auto rho = hamil.density_matrix(point.first, point.second);
        double time = point.first / decay_rate;
        write_state_info(rho_out, time, rho, hamil);
        write_kdist(kdistout, time, rho, hamil);
    }
    rho_out.close();
    kdistout.close();

    // Output just the final k distribution for convenience
    kdistfinalout << kdist_header << std::endl;
    auto rho = hamil.density_matrix(
        rho_c_solution.back().first, rho_c_solution.back().second);
    write_kdist(kdistfinalout, rho_c_solution.back().first / decay_rate,
        rho, hamil);
}

void write_state_info(std::ofstream& outfile, double t,
    const std::vector<std::complex<double>>& rho, const HMotion& hamil) {
    outfile << t
        << " " << std::real(hamil.partialtr_k(rho, 0))
        << " " << std::real(hamil.partialtr_k(rho, 1))
        << " " << std::real(hamil.partialtr_k(rho, 2))
        << " " << std::real(hamil.totaltr(rho))
        << " " << std::real(hamil.purity(rho));
    double krms = 0;    // RMS k value
    for(int k = 0; k <= hamil.kmax; ++k) {
        auto ktr = hamil.partialtr_n(rho, k);
        if(k != 0) {
            ktr += hamil.partialtr_n(rho, -k);
        }
        outfile << " " << std::real(ktr);
        krms += std::real(ktr) * k*k;
    }
    outfile << " " << sqrt(krms);
    outfile << std::endl;
}
void write_kdist(std::ofstream& outfile, double t,
    const std::vector<std::complex<double>>& rho, const HMotion& hamil) {
    for(int k = -hamil.kmax; k <= hamil.kmax; ++k) {
        outfile << t
            << " " << k << " " << std::real(hamil.partialtr_n(rho, k))
            << " " << std::real(rho[hamil.subidx(0, k, 0, k)])
            << " " << std::real(rho[hamil.subidx(1, k, 1, k)])
            << " " << std::real(rho[hamil.subidx(2, k, 2, k)])
            << std::endl;
    }
}