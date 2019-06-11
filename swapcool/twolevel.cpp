// The density matrix is stored in row major format.
#include "twolevel.hpp"

const std::string DEFAULT_CFG_FILE = "params.cfg";
const std::string OUTFILEBASE = "output/rho_twolevel.out";
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
    HSawtoothRotWave deriv(cfg_file);

    std::ofstream cyclesout("cycles.out");
    for(double gt = 0; gt < duration_by_decay; gt += 1e-3) {
        cyclesout << gt << " " << deriv.rabi_softswitch(gt) << " "
        << deriv.detun_per_decay(gt) << " " << deriv.cumulative_phase(gt)
        << std::endl;
    }
    cyclesout.close();

    // Set the initial condition to be all in the ground state
    std::vector<std::complex<double>> rho_c0{1, 0, 0, 0};

    // Solve the system in natural units with an adaptive RK method
    auto rho_c_solution = timestepping::odesolve(deriv, rho_c0,
        duration_by_decay, timestepping::AdaptiveRK(tol));

    // Form filename
    std::ostringstream oftag_ss;
    oftag_ss << std::setprecision(OUTFILENAME_PRECISION)
        << "A" << deriv.detun_amp_per_decay
        << "_f" << deriv.detun_freq_per_decay
        << "_Omega" << deriv.rabi_freq_per_decay;
    std::string ofname = OUTFILEBASE;
    ofname.insert(ofname.rfind("."), "_" + oftag_ss.str());

    // Write the solution in SI units, putting back in the rotating wave
    // oscillation
    // Column order is t, |rho11|, |rho22|
    std::ofstream outfile(ofname);
    // Write table header
    outfile << "t |rho11| |rho22|" << std::endl;
    for(auto point: rho_c_solution) {
        // Convert rho_c (rotating wave coefficients) to rho
        auto rho = deriv.density_matrix(point.first, point.second);
        outfile << point.first / decay_rate // Convert from Gamma*t to just t
            << " " << std::abs(rho[0]) << " " << std::abs(rho[3]) << std::endl;
    }
    outfile.close();
}
