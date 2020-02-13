[Paper](https://doi.org/10.1103/PhysRevA.101.013422) ([e-Print](https://arxiv.org/abs/1909.11698))

A collection of simulations related to different laser cooling methods.

# Usage
Run `make` in the top-level directory. The generated executables are stored in the `bin/` directory.

# Simulations
- `optical_molasses` (source code in `src/optmol/`)  contains semiclassical Monte Carlo simulations for the standard "optical molasses" laser cooling.
- `swapint` and `swapmotion` (source code in `src/swapcool/`) contains density matrix simulations for Sawtooth-Wave Adiabatic Passage (SWAP) cooling, described by [Bartolotta et. al., Physical Review A 98, 023404 (2018)](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.98.023404). `swapint` is a simulation of just internal states, while `swapmotion` accounts for both internal and momentum states.

# Notes on some other directories
- `config/` holds default configuration files for the simulations.
- `doc/` holds explanations of the physics of the simulated cooling methods, as well as code documentation.
- `scripts/plotting/` contains simple Python scripts for quick and dirty plotting.
- `vendor/pcg-cpp-0.98/` contains the PCG RNG, which provides faster random number generation than the C++ <random> library.

# Acknowledgements

Research is done under Professor Francis Robicheaux for the 2019 Research Experience for Undergraduates (REU) program at Purdue University.
