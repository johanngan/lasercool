// SWAP laser setup only considering internal states, without particle motion
// The density matrix is stored in row major format.
#include "swapint.hpp"

const std::string DEFAULT_CFG_FILE = "params.cfg";
const std::string OUTFILEBASE = "output/rho_swapint.out";
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
    HInt hamil(cfg_file);

    std::ofstream cyclesout("cycles.out");
    for(double gt = 0; gt < duration_by_decay; gt += 1e-3) {
        cyclesout << gt << " " << hamil.rabi_softswitch(gt) << " "
        << hamil.detun_per_decay(gt) << " " << hamil.cumulative_phase(gt)
        << std::endl;
    }
    cyclesout.close();

    // Set the initial condition to be all in the ground state
    std::vector<std::complex<double>> rho_c0{
        0, 0, 0,
        0, 1, 0,
        0, 0, 0,
    };

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

    // Write the solution in SI units, putting back in the rotating wave
    // oscillation
    // Column order is t, |rho11|, |rho22|
    std::ofstream outfile(ofname);
    // Write table header
    outfile << "t |rho11| |rho22| |rho33|" << std::endl;
    for(auto point: rho_c_solution) {
        // Convert rho_c (rotating wave coefficients) to rho
        auto rho = hamil.density_matrix(point.first, point.second);
        outfile << point.first / decay_rate // Convert from Gamma*t to just t
            << " " << std::abs(rho[hamil.subidx(0,0)])
            << " " << std::abs(rho[hamil.subidx(1,1)])
            << " " << std::abs(rho[hamil.subidx(2,2)])
            << std::endl;
    }
    outfile.close();
}
