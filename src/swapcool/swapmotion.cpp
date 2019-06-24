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

    std::cout << "In units of decay rate when applicable:" << std::endl
        << "    Decay rate: " << decay_rate << std::endl
        << "    Decay: " << (hamil.enable_decay ? "on" : "off") << std::endl
        << "    Branching ratio: " << hamil.branching_ratio << std::endl
        << "    Delta amplitude: " << hamil.detun_amp_per_decay << std::endl
        << "    Sawtooth frequency: " << hamil.detun_freq_per_decay << std::endl
        << "    Rabi frequency: " << hamil.rabi_freq_per_decay << std::endl
        << "    Recoil frequency: " << hamil.recoil_freq_per_decay << std::endl
        << "    Initial momentum state: " << init_k << std::endl
        << "    Maximum momentum state: " << hamil.kmax << std::endl
        << "    Stepper tolerance: " << tol << std::endl;

    // Form output files
    std::ostringstream oftag_ss;
    oftag_ss << std::setprecision(OUTFILENAME_PRECISION)
        << "A" << hamil.detun_amp_per_decay
        << "_f" << hamil.detun_freq_per_decay
        << "_Omega" << hamil.rabi_freq_per_decay
        << "_recoil" << hamil.recoil_freq_per_decay
        << "_" << (hamil.enable_decay ? "" : "no") << "decay"
        << "_B" << hamil.branching_ratio
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

    // Write table headers
    rho_out << "t |rho11| |rho22| |rho33| tr(rho) tr(rho^2)";
    for(int k = 0; k <= hamil.kmax; ++k) {
        rho_out << " |k" << k << "|";
    }
    rho_out << " |k_rms|";
    rho_out << std::endl;
    std::string kdist_header = "t k P(k) P(n = 0, k), P(n = 1, k), P(n = 2, k)";
    kdistout << kdist_header << std::endl;
    kdistfinalout << kdist_header << std::endl;

    // Solve the system
    // Figure out how many cycles to run.
    double nfullcycles_double;
    double cycle_remain = modf(
        hamil.detun_freq_per_decay*duration_by_decay, &nfullcycles_double);
    int nfullcycles = static_cast<int>(nfullcycles_double);
    bool has_partial_cycle = (cycle_remain != 0);

    std::vector<std::complex<double>> rho_c(hamil.nmat());
    // Initialize all in one k-state
    rho_c[hamil.subidx(1, init_k, 1, init_k)] = 1;
    // For holding the time of the popped final entry of the solution,
    // to be used after loop termination
    double solution_endgt = 0;
    // Solve cycle-by-cycle. Add an extra iteration if a partial cycle is
    // necessary
    for(int cycle = 0; cycle < nfullcycles + has_partial_cycle; ++cycle) {
        // Determine the final local cycle time to solve until
        double endtime = std::min(
            duration_by_decay, (cycle+1)/hamil.detun_freq_per_decay)
            - cycle/hamil.detun_freq_per_decay;
        
        // Prepare the density matrix for a new cycle
        initialize_cycle(rho_c, hamil);
        // Solve a full/partial system cycle in natural units with adaptive RK
        auto rho_c_solution = timestepping::odesolve(hamil, rho_c,
            endtime, timestepping::AdaptiveRK(tol));

        // Save the final rho_c for the next cycle        
        std::tie(solution_endgt, rho_c) = rho_c_solution.back();
        solution_endgt += cycle/hamil.detun_freq_per_decay;
        // Don't write the final state to file, since it'll be modified and
        // included in the next iteration, or written after loop exit
        rho_c_solution.pop_back();
        
        // Write the solution to file
        for(auto point: rho_c_solution) {
            // Get the actual, global time
            double gt = point.first + cycle/hamil.detun_freq_per_decay;
            double time = gt / decay_rate;

            auto rho = hamil.density_matrix(gt, point.second);
            write_state_info(rho_out, time, rho, hamil);
            write_kdist(kdistout, time, rho, hamil);
        }

        std::cout << "Completed cycle " << cycle << std::endl;
    }
    // Write the final state to file
    double solution_endtime = solution_endgt / decay_rate;
    auto rhofinal = hamil.density_matrix(solution_endtime, rho_c);
    write_state_info(rho_out, solution_endtime, rhofinal, hamil);
    write_kdist(kdistout, solution_endtime, rhofinal, hamil);
    rho_out.close();
    kdistout.close();

    // Output just the final k distribution to a separate file for convenience
    write_kdist(kdistfinalout, solution_endtime, rhofinal, hamil);
    kdistfinalout.close();
}

void initialize_cycle(std::vector<std::complex<double>>& rho,
    const HMotion& hamil) {
    for(int kl = -hamil.kmax; kl <= hamil.kmax; ++kl) {
        for(int kr = -hamil.kmax; kr <= hamil.kmax; ++kr) {
            // Excited state population and intra-excited-state coherences
            // distribute between the lower energy states
            rho[hamil.subidx(0, kl, 0, kr)] += (1 - hamil.branching_ratio)
                * rho[hamil.subidx(2, kl, 2, kr)];
            rho[hamil.subidx(1, kl, 1, kr)] +=
                hamil.stationary_decay_prob*hamil.branching_ratio
                * rho[hamil.subidx(2, kl, 2, kr)];
            if(kl - 1 >= -hamil.kmax && kr - 1 >= -hamil.kmax) {
                rho[hamil.subidx(1, kl-1, 1, kr-1)] +=
                    (1-hamil.stationary_decay_prob)/2*hamil.branching_ratio
                    * rho[hamil.subidx(2, kl, 2, kr)];
            }
            if(kl + 1 <= hamil.kmax && kr + 1 <= hamil.kmax) {
                rho[hamil.subidx(1, kl+1, 1, kr+1)] +=
                    (1-hamil.stationary_decay_prob)/2*hamil.branching_ratio
                    * rho[hamil.subidx(2, kl, 2, kr)];
            }

            // Excited state and excited-state coherences decay to 0
            rho[hamil.subidx(2, kl, 2, kr)] = 0;
            for(unsigned n = 0; n < 2; ++n) {
                rho[hamil.subidx(n, kl, 2, kr)] = 0;
                rho[hamil.subidx(2, kl, n, kr)] = 0;
            }
        }
    }
    return;
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