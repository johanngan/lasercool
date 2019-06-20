// SWAP laser setup considering both internal states and momentum states.
// The density matrix is stored in row major format, enumerated
// as |nl, kl><nr, kr|
#include "swapmotion.hpp"

const std::string DEFAULT_CFG_FILE = "config/params_swapcool.cfg";
const std::string OUTFILEBASE = "output/swapcool/rho_swapmotion.out";
const unsigned OUTFILENAME_PRECISION = 3;

int main(int argc, char** argv) {
    // Read in a possible config file
    std::string cfg_file = DEFAULT_CFG_FILE;
    if(argc > 1) {
        cfg_file = std::string(argv[1]);
    }

    double decay_rate, duration_by_decay, tol;
    load_params(cfg_file,
        {
            {"spontaneous_decay_rate", &decay_rate},
            {"duration", &duration_by_decay},
            {"tolerance", &tol}
        }
    );

    // Form the derivative operator, in natural units
    // d(rho)/d(Gamma*t)
    HMotion hamil(cfg_file);

    std::vector<std::complex<double>> rho_c0(hamil.nmat());
    int k0 = 0;
    rho_c0[hamil.subidx(1, k0, 1, k0)] = 1;

    // Solve the system in natural units with an adaptive RK method
    auto rho_c_solution = timestepping::odesolve(hamil, rho_c0,
        duration_by_decay, timestepping::AdaptiveRK(tol));

    // Form filename
    std::ostringstream oftag_ss;
    oftag_ss << std::setprecision(OUTFILENAME_PRECISION)
        << "A" << hamil.detun_amp_per_decay
        << "_f" << hamil.detun_freq_per_decay
        << "_Omega" << hamil.rabi_freq_per_decay;
    std::string ofname = OUTFILEBASE;
    ofname.insert(ofname.rfind("."), "_" + oftag_ss.str());

    // Write the solution in SI units
    // Column order is t, |rho11|, |rho22|, |rho33|
    std::ofstream outfile(ofname);
    // Write table header
    outfile << "t |rho11| |rho22| |rho33| tr(rho) tr(rho^2)";
    for(int k = 0; k <= hamil.kmax; ++k) {
        outfile << " |k" << k << "|";
    }
    outfile << " |k_rms|";
    outfile << std::endl;
    for(auto point: rho_c_solution) {
        auto rho = hamil.density_matrix(point.first, point.second);
        outfile << point.first / decay_rate // Convert from Gamma*t to just t
            << " " << std::abs(hamil.partialtr_k(rho, 0))
            << " " << std::abs(hamil.partialtr_k(rho, 1))
            << " " << std::abs(hamil.partialtr_k(rho, 2))
            << " " << std::abs(hamil.totaltr(rho))
            << " " << std::abs(hamil.purity(rho));
            double krms = 0;    // RMS k value
            for(int k = 0; k <= hamil.kmax; ++k) {
                auto ktr = hamil.partialtr_n(rho, k);
                if(k != 0) {
                    ktr += hamil.partialtr_n(rho, -k);
                }
                outfile << " " << std::abs(ktr);
                krms += std::abs(ktr) * k*k;
            }
            outfile << " " << sqrt(krms);
            outfile << std::endl;
    }
    outfile.close();

    std::ofstream khist("output/swapcool/khist.out");
    auto rhof = hamil.density_matrix(
        rho_c_solution.back().first, rho_c_solution.back().second);
    for(int k = -hamil.kmax; k <= hamil.kmax; ++k) {
        auto ktr = hamil.partialtr_n(rhof, k);
        khist << k << " " << std::abs(ktr)
            << " " << k*k << " " << std::abs(ktr)*std::abs(ktr) << std::endl;
    }
    khist.close();
}