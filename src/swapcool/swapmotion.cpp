// SWAP laser setup considering both internal states and momentum states.
// The density matrix is stored in row major format, enumerated
// as |nl, kl><nr, kr|
#include "swapmotion.hpp"

const std::string DEFAULT_CFG_FILE = "config/params_swapcool.cfg";
const std::string OUTPUT_DIR = "output/swapcool/swapmotion";
const std::string RHO_OUTFILEBASE = "rho.out";
const std::string KDIST_OUTFILEBASE = "kdist.out";
const std::string KDIST_FINAL_OUTFILEBASE = "kdist_final.out";
// Approximate number of solution points to output per sawtooth cycle.
// Only approximate because adaptive time steps make it hard to divide things
// exactly
const double APPROX_OUTPUT_PTS_PER_CYCLE = 100;
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

    // Print out parameters
    std::cout << "In units of decay rate when applicable:" << std::endl
        << "    Decay rate: " << decay_rate << std::endl
        << "    Decay: " << (hamil.enable_decay ? "on" : "off") << std::endl
        << "    Branching ratio: " << hamil.branching_ratio << std::endl
        << "    Delta amplitude: " << hamil.detun_amp_per_decay << std::endl
        << "    Sawtooth frequency: " << hamil.detun_freq_per_decay << std::endl
        << "    Rabi frequency: " << hamil.rabi_freq_per_decay << std::endl
        << "    Recoil frequency: " << hamil.recoil_freq_per_decay << std::endl
        << "    Initial momentum state: " << init_k << std::endl
        << "    Momentum state range: ["
            << hamil.handler.kmin << ", " << hamil.handler.kmax
        << "]" << std::endl
        << "    Duration: " << duration_by_decay << " ("
        << hamil.detun_freq_per_decay*duration_by_decay << " cycles)"
        << std::endl
        << "    Stepper tolerance: " << tol << std::endl
        << std::endl;
    
    // Calculate quality metrics
    double dopshift = hamil.recoil_freq_per_decay*init_k;
    double rampsize = hamil.detun_amp_per_decay / (4*dopshift);
    double qfactor = 0.5*hamil.detun_amp_per_decay*hamil.detun_freq_per_decay
        / (dopshift - hamil.recoil_freq_per_decay + hamil.rabi_freq_per_decay);
    double adiabaticity = hamil.rabi_freq_per_decay*hamil.rabi_freq_per_decay
        / (2*hamil.detun_amp_per_decay*hamil.detun_freq_per_decay);
    double splitting = 2*(dopshift - hamil.recoil_freq_per_decay)
        / hamil.rabi_freq_per_decay;
    
    std::cout << "Quality metrics (* mildly low, ** very low):" << std::endl
        << "    " << evaluate_quality_metric(rampsize)
        << "Ramp size: " << rampsize << std::endl
        << "    " << evaluate_quality_metric(qfactor)
        << "Q factor: " << qfactor << std::endl
        << "    " << evaluate_quality_metric(adiabaticity)
        << "Adiabaticity: " << adiabaticity << std::endl
        << "    " << evaluate_quality_metric(splitting)
        << "Doppler splitting: " << splitting << std::endl
        << std::endl;

    // Approximate gamma*dt between output points
    double output_gdt = 1. /
        (APPROX_OUTPUT_PTS_PER_CYCLE * hamil.detun_freq_per_decay);

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
    for(int k = (hamil.handler.kmax*hamil.handler.kmin <= 0 ?
            0 : std::min(std::abs(hamil.handler.kmax),
                         std::abs(hamil.handler.kmin)));
        k <= std::max(std::abs(hamil.handler.kmax),
                      std::abs(hamil.handler.kmin));
        ++k) {
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

    std::vector<std::complex<double>> rho_c(hamil.handler.idxmap.size());
    // Initialize all in one k-state
    hamil.handler.at(rho_c, 1, init_k, 1, init_k) = 1;
    // For holding the time of the popped final entry of the solution,
    // to be used after loop termination
    double solution_endgt = 0;
    // Solve cycle-by-cycle. Add an extra iteration if a partial cycle is
    // necessary

    /// TIMING
    auto start = std::chrono::system_clock::now();
    ///

    for(int cycle = 0; cycle < nfullcycles + has_partial_cycle; ++cycle) {
        std::cout << "\rProgress: running cycle " << cycle + 1
            << "/" << nfullcycles + has_partial_cycle << std::flush;

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
        int cur_steps = -1; // Effective number of output steps taken so far
        for(auto point: rho_c_solution) {
            // Get the actual, global time
            double gt = point.first + cycle/hamil.detun_freq_per_decay;
            double time = gt / decay_rate;

            // Get the effective number of output time steps taken so far
            int cur_steps_new = static_cast<int>(gt / output_gdt);
            // Only record output if time has advanced by at least the minimum
            // specified time between outputs
            if(cur_steps_new > cur_steps) {
                // Record the new number of output time steps taken
                cur_steps = cur_steps_new;

                auto rho = hamil.density_matrix(gt, point.second);
                write_state_info(rho_out, time, rho, hamil.handler);
                write_kdist(kdistout, time, rho, hamil.handler);
            }
        }
    }
    std::cout << std::endl;

    ///
    std::chrono::duration<double> total_seconds =
        std::chrono::system_clock::now() - start;
    std::cout << "Simulation time: " << total_seconds.count() << " s"
        << std::endl;
    ///

    // Write the final state to file
    double solution_endtime = solution_endgt / decay_rate;
    auto rhofinal = hamil.density_matrix(solution_endtime, rho_c);
    write_state_info(rho_out, solution_endtime, rhofinal, hamil.handler);
    write_kdist(kdistout, solution_endtime, rhofinal, hamil.handler);
    rho_out.close();
    kdistout.close();

    // Output just the final k distribution to a separate file for convenience
    write_kdist(kdistfinalout, solution_endtime, rhofinal, hamil.handler);
    kdistfinalout.close();
}

std::string evaluate_quality_metric(
    double metric,
    double low_thresh,
    double very_low_thresh,
    std::string okay_str,
    std::string low_str,
    std::string very_low_str) {
    if(metric < very_low_thresh) {
        return very_low_str;
    } else if(metric < low_thresh) {
        return low_str;
    }
    return okay_str;
}

void initialize_cycle(std::vector<std::complex<double>>& rho,
    const HMotion& hamil) {
    for(int kl = hamil.handler.kmin; kl <= hamil.handler.kmax; ++kl) {
        for(int kr = hamil.handler.kmin; kr <= hamil.handler.kmax; ++kr) {
            // Excited state population and intra-excited-state coherences
            // distribute between the lower energy states
            if(hamil.handler.has(0, kl, 0, kr)) {
                hamil.handler.at(rho, 0, kl, 0, kr) += (1 - hamil.branching_ratio)
                    * hamil.handler.ele(rho, 2, kl, 2, kr);
            }
            if(hamil.handler.has(1, kl, 1, kr)) {
                hamil.handler.at(rho, 1, kl, 1, kr) +=
                    hamil.stationary_decay_prob*hamil.branching_ratio
                    * hamil.handler.ele(rho, 2, kl, 2, kr);
            }
            if(kl - 1 >= hamil.handler.kmin && kr - 1 >= hamil.handler.kmin
                && hamil.handler.has(1, kl-1, 1, kr-1)) {
                hamil.handler.at(rho, 1, kl-1, 1, kr-1) +=
                    (1-hamil.stationary_decay_prob)/2*hamil.branching_ratio
                    * hamil.handler.ele(rho, 2, kl, 2, kr);
            }
            if(kl + 1 <= hamil.handler.kmax && kr + 1 <= hamil.handler.kmax
                && hamil.handler.has(1, kl+1, 1, kr+1)) {
                hamil.handler.at(rho, 1, kl+1, 1, kr+1) +=
                    (1-hamil.stationary_decay_prob)/2*hamil.branching_ratio
                    * hamil.handler.ele(rho, 2, kl, 2, kr);
            }

            // Excited state and excited-state coherences decay to 0
            if(hamil.handler.has(2, kl, 2, kr)) {
                hamil.handler.at(rho, 2, kl, 2, kr) = 0;
            }
            if(hamil.handler.has(1, kl, 2, kr)) {
                hamil.handler.at(rho, 1, kl, 2, kr) = 0;
            }
            if(hamil.handler.has(2, kl, 1, kr)) {
                hamil.handler.at(rho, 2, kl, 1, kr) = 0;
            }
        }
    }
    return;
}

void write_state_info(std::ofstream& outfile, double t,
    const std::vector<std::complex<double>>& rho, const DensMatHandler& handler) {
    outfile << t;
    for(unsigned n = 0; n < handler.nint; ++n) {
        outfile << " " << std::real(handler.partialtr_k(rho, n));
    }
    outfile << " " << std::real(handler.totaltr(rho))
        << " " << std::real(handler.purity(rho));
    double krms = 0;    // RMS k value
    for(int k = handler.kmin; k <= handler.kmax; ++k) {
        krms += std::real(handler.partialtr_n(rho, k))*k*k;
    }
    outfile << " " << sqrt(krms);
    outfile << std::endl;
}
void write_kdist(std::ofstream& outfile, double t,
    const std::vector<std::complex<double>>& rho, const DensMatHandler& handler) {
    for(int k = handler.kmin; k <= handler.kmax; ++k) {
        outfile << t << " " << k
            << " " << std::real(handler.partialtr_n(rho, k));
        for(unsigned n = 0; n < handler.nint; ++n) {
            outfile << " " << std::real(handler.ele(rho, n, k, n, k));
        }
        outfile << std::endl;
    }
}