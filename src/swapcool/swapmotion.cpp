// SWAP laser setup considering both internal states and momentum states.
// The density matrix is stored in row major format, enumerated
// as |nl, kl><nr, kr|
#include "swapmotion.hpp"

const std::string DEFAULT_CFG_FILE = "config/params_swapcool.cfg";
const std::string DEFAULT_OUTPUT_DIR = "output/swapcool/swapmotion";
const std::string RHO_OUTFILEBASE = "rho.out";
const std::string KDIST_OUTFILEBASE = "kdist.out";
const std::string KDIST_FINAL_OUTFILEBASE = "kdist_final.out";
// Approximate number of solution points to output per sawtooth cycle.
// Only approximate because adaptive time steps make it hard to divide things
// exactly
const double APPROX_OUTPUT_PTS_PER_CYCLE = 100;
const unsigned OUTFILENAME_PRECISION = 3;

template<typename T>
inline T sqr(T x) {return x*x;}

int main(int argc, char** argv) {
    // Parse the program name to find the project root directory
    std::string progdir, progname;
    std::tie(progname, progdir) = fileparts(argv[0]);
    // The program binary will be in project/bin, assuming no symlinks
    std::string projrootdir = progdir + "/..";

    if(argc > 4) {
        std::cout << "Usage: " << progname
            << " [<output directory>] [<config file>] [--batch-mode]"
            << std::endl;
        return 1;
    }
    // Read in a possible output directory
    std::string output_dir = fullfile(DEFAULT_OUTPUT_DIR, projrootdir);
    if(argc > 1) {
        output_dir = std::string(argv[1]);
    }
    // Read in a possible config file
    std::string cfg_file = fullfile(DEFAULT_CFG_FILE, projrootdir);
    if(argc > 2) {
        cfg_file = std::string(argv[2]);
    }
    // In "batch mode", don't output any info to the console
    bool batchmode = false;
    if(argc > 3) {
        std::string mode(argv[3]);
        if(mode == "-b" || mode == "--batch-mode") {
            batchmode = true;
        } else {
            std::cout << "Invalid argument: "
                "\"-b\" or \"--batch-mode\" for batch mode, "
                "otherwise omit third argument." << std::endl;
            return 1;
        }
    }

    double duration_by_decay, tol, init_temp, init_k;
    load_params(cfg_file,
        {
            {"duration", &duration_by_decay},
            {"tolerance", &tol},
            {"initial_temperature", &init_temp},
            {"initial_momentum", &init_k}
        }
    );
    bool is_thermal = true;
    if(!std::isnan(init_k)) {
        // Override temperature and start from a fixed k
        is_thermal = false;
    }

    // Form the derivative operator, in natural units
    // d(rho)/d(Gamma*t)
    HMotion hamil(cfg_file);

    // Initialize state
    std::vector<std::complex<double>> rho_c;
    if(is_thermal) {
        // Initialize to thermal distribution
        rho_c = thermal_state(init_temp, hamil);
    } else {
        // Initialize all in one k-state
        int kidx = hamil.handler.closest_kidx(init_k);
        rho_c.resize(hamil.handler.idxmap.size());
        hamil.handler.at(rho_c, 1, kidx, 1, kidx) = 1;
    }

    // Print out stuff if not in batch mode
    if(!batchmode) {
        print_system_info(rho_c, hamil, init_temp, init_k, is_thermal,
            duration_by_decay, tol);
    }

    // Form output files
    std::ostringstream oftag_ss;
    oftag_ss << std::setprecision(OUTFILENAME_PRECISION)
        << "A" << hamil.detun_amp_per_decay
        << "_f" << hamil.detun_freq_per_decay
        << "_Omega" << hamil.rabi_freq_per_decay
        << "_recoil" << hamil.recoil_freq_per_decay
        << "_" << (hamil.enable_decay ? "" : "no") << "decay"
        << "_B" << hamil.branching_ratio;
    if(is_thermal) {
        oftag_ss << "_T" << init_temp;
    } else {
        oftag_ss << "_k" << init_k;
    }

    std::ofstream rho_out(fullfile(tag_filename(
        RHO_OUTFILEBASE, oftag_ss.str()),
        output_dir
    ));
    std::ofstream kdistout(fullfile(tag_filename(
        KDIST_OUTFILEBASE, oftag_ss.str()),
        output_dir
    ));
    std::ofstream kdistfinalout(fullfile(tag_filename(
        KDIST_FINAL_OUTFILEBASE, oftag_ss.str()),
        output_dir
    ));

    // Write table headers
    rho_out << "t |rho11| |rho22| |rho33| tr(rho) tr(rho^2)"
        << " |k_rms| |k_rms(unleaked)|" << std::endl;
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

    // Approximate gamma*dt between output points
    double output_gdt = 1. /
        (APPROX_OUTPUT_PTS_PER_CYCLE * hamil.detun_freq_per_decay);
    // For holding the time of the popped final entry of the solution,
    // to be used after loop termination
    double solution_endgt = 0;
    // Solve cycle-by-cycle. Add an extra iteration if a partial cycle is
    // necessary

    /// TIMING
    auto start = std::chrono::system_clock::now();
    ///

    for(int cycle = 0; cycle < nfullcycles + has_partial_cycle; ++cycle) {
        if(!batchmode) {
            std::cout << "\rProgress: running cycle " << cycle + 1
                << "/" << nfullcycles + has_partial_cycle << std::flush;
        }

        // Determine the final local cycle time to solve until
        double endtime = std::min(
            duration_by_decay, (cycle+1)/hamil.detun_freq_per_decay)
            - cycle/hamil.detun_freq_per_decay;
        
        // Prepare the density matrix for a new cycle
        hamil.initialize_cycle(rho_c);
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
            double time = gt / hamil.decay_rate;

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
    if(!batchmode) {
        std::cout << std::endl;

        ///
        std::chrono::duration<double> total_seconds =
            std::chrono::system_clock::now() - start;
            std::cout << "Simulation time: " << total_seconds.count() << " s"
            << std::endl;
        ///
    }

    // Write the final state to file
    double solution_endtime = solution_endgt / hamil.decay_rate;
    auto rhofinal = hamil.density_matrix(solution_endtime, rho_c);
    write_state_info(rho_out, solution_endtime, rhofinal, hamil.handler);
    write_kdist(kdistout, solution_endtime, rhofinal, hamil.handler);
    rho_out.close();
    kdistout.close();

    // Output just the final k distribution to a separate file for convenience
    write_kdist(kdistfinalout, solution_endtime, rhofinal, hamil.handler);
    kdistfinalout.close();
}

std::vector<std::complex<double>> thermal_state(double temp,
    const HMotion& hamil) {
    std::vector<std::complex<double>> rho(hamil.handler.idxmap.size());
    double partition_fn = 0;
    for(int kidx = hamil.handler.kmin; kidx <= hamil.handler.kmax; ++kidx) {
        double boltz_weight = std::exp(-fundamental_constants::HBAR
            *hamil.recoil_freq_per_decay*hamil.decay_rate
            * sqr(hamil.handler.kval(kidx))
            / (fundamental_constants::K_BOLTZMANN*temp));
        partition_fn += boltz_weight;
        hamil.handler.at(rho, 1, kidx, 1, kidx) = boltz_weight;
    }
    // Normalize by partition function
    for(int kidx = hamil.handler.kmin; kidx <= hamil.handler.kmax; ++kidx) {
        hamil.handler.at(rho, 1, kidx, 1, kidx) /= partition_fn;
    }
    return rho;
}

void print_system_info(const std::vector<std::complex<double>>& rho,
    const HMotion& hamil, double init_temp, double init_k, bool is_thermal,
    double duration_by_decay, double tol) {
    // Parameters
    std::cout << "In units of decay rate when applicable:" << std::endl
        << "    Decay rate: " << hamil.decay_rate << std::endl
        << "    Decay: " << (hamil.enable_decay ? "on" : "off")
        << std::endl
        << "    Branching ratio: " << hamil.branching_ratio << std::endl
        << "    Delta amplitude: " << hamil.detun_amp_per_decay
        << std::endl
        << "    Sawtooth frequency: " << hamil.detun_freq_per_decay
        << std::endl
        << "    Rabi frequency: " << hamil.rabi_freq_per_decay << std::endl
        << "    Recoil frequency: " << hamil.recoil_freq_per_decay
        << std::endl;

    if(is_thermal) {
        std::cout << "    Initial temperature: " << init_temp << " K"
            << std::endl;
    } else {
        std::cout << "    Initial momentum state: " << init_k << std::endl;
    }
    std::cout << "    Momentum state range: ["
        << hamil.handler.kval(hamil.handler.kmin) << ", "
        << hamil.handler.kval(hamil.handler.kmax) << "]" << std::endl
        << "    Momentum state subdivisions: " << hamil.handler.ksubdivs
        << std::endl
        << "    Duration: " << duration_by_decay << " ("
        << hamil.detun_freq_per_decay*duration_by_decay << " cycles)"
        << std::endl
        << "    Stepper tolerance: " << tol << std::endl
        << std::endl;

    // Quality metrics
    double dopshift = hamil.recoil_freq_per_decay
        * calc_krms(rho, hamil.handler);
    double rampsize = hamil.detun_amp_per_decay / (4*dopshift);
    double qfactor = hamil.detun_amp_per_decay*hamil.detun_freq_per_decay
        / (2*(dopshift - hamil.recoil_freq_per_decay) + hamil.rabi_freq_per_decay);
    double adiabaticity = hamil.rabi_freq_per_decay*hamil.rabi_freq_per_decay
        / (2*hamil.detun_amp_per_decay*hamil.detun_freq_per_decay);
    double splitting = 2*(dopshift - hamil.recoil_freq_per_decay)
        / hamil.rabi_freq_per_decay;

    std::cout << "Quality metrics (* mildly low, ** very low):"
        << std::endl
        << "    " << evaluate_quality_metric(rampsize)
        << "Ramp size: " << rampsize << std::endl
        << "    " << evaluate_quality_metric(qfactor)
        << "Q factor: " << qfactor << std::endl
        << "    " << evaluate_quality_metric(adiabaticity)
        << "Adiabaticity: " << adiabaticity << std::endl
        << "    " << evaluate_quality_metric(splitting)
        << "Doppler splitting: " << splitting << std::endl
        << std::endl;
}

double calc_krms(const std::vector<std::complex<double>>& rho,
    const DensMatHandler& handler) {
    double krms = 0;
    for(int kidx = handler.kmin; kidx <= handler.kmax; ++kidx) {
        krms += std::real(handler.partialtr_n(rho, kidx))*sqr(handler.kval(kidx));
    }
    return sqrt(krms);
}

double calc_krms_unleaked(const std::vector<std::complex<double>>& rho,
    const DensMatHandler& handler) {
    // For renormalization
    double unleaked_prob = std::real(
        handler.partialtr_k(rho, 1) + handler.partialtr_k(rho, 2));

    double krms = 0;
    for(int kidx = handler.kmin; kidx <= handler.kmax; ++kidx) {
        double prob = std::real(handler.at(rho, 1, kidx, 1, kidx) \
            + handler.at(rho, 2, kidx, 2, kidx)) / unleaked_prob;
        krms += prob * sqr(handler.kval(kidx));
    }
    return sqrt(krms);
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

void write_state_info(std::ofstream& outfile, double t,
    const std::vector<std::complex<double>>& rho, const DensMatHandler& handler) {
    outfile << t;
    for(unsigned n = 0; n < handler.nint; ++n) {
        outfile << " " << std::real(handler.partialtr_k(rho, n));
    }
    outfile << " " << std::real(handler.totaltr(rho))
        << " " << std::real(handler.purity(rho))
        << " " << calc_krms(rho, handler)
        << " " << calc_krms_unleaked(rho, handler);
    outfile << std::endl;
}
void write_kdist(std::ofstream& outfile, double t,
    const std::vector<std::complex<double>>& rho, const DensMatHandler& handler) {
    for(int kidx = handler.kmin; kidx <= handler.kmax; ++kidx) {
        outfile << t << " " << handler.kval(kidx)
            << " " << std::real(handler.partialtr_n(rho, kidx));
        for(unsigned n = 0; n < handler.nint; ++n) {
            outfile << " " << std::real(handler.ele(rho, n, kidx, n, kidx));
        }
        outfile << std::endl;
    }
}