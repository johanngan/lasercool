A collection of simulations related to different laser cooling methods.

# Simulations
- `optmol/` contains semiclassical Monte Carlo simulations for the standard "optical molasses" laser cooling.
- `swapcool/` contains density matrix simulations for Sawtooth-Wave Adiabatic Passage (SWAP) cooling, described by [Bartolotta et. al., Physical Review A 98, 023404 (2018)](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.98.023404).

# Other directories
- `lib/` contains reusable utility code for the different simulations.
- `plotting/` contains simple Python scripts for quick and dirty plotting
- `vendor/pcg-cpp-0.98/` contains the PCG RNG, which provides faster random number generation than the C++ <random> library.

# Acknowledgements

Research is done under Professor Francis Robicheaux for the 2019 Research Experience for Undergraduates (REU) program at Purdue University.
