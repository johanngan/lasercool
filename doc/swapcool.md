Density matrix simulation of the "sawtooth-wave adiabatic passage" laser cooling technique.

# General explanation
TODO

# Simulation details
TODO

## Detuning sawtooth oscillation
TODO

## Internal state simulation
TODO

## Motional and internal state simulation
TODO

# Usage
Run `make swapcool` in the top-level directory, set the parameters in `/config/params_swapcool.cfg`, then run `/bin/swapint` or `/bin/swapmotion`. Optionally give the path to the directory to write output to, and the path to the configuration file to use. For `swapmotion`, a final optional parameter, `--batch-mode` (or `-b` for short) can be specified to enable batch mode, which suppresses all console output.

## OpenMP Capability
If OpenMP is available on your machine, enable it by adding the appropriate compiler/linker flags when running make. I.e. compile swapcool with `make swapcool CFLAGS=-openmp FLAGS=-fopenmp`.

Set the number of threads with the environment variable `OMP_NUM_THREADS`. E.g. specify 4 threads by running `export OMP_NUM_THREADS=4`.

## Lab parameters
`params_swapcool.cfg` contains different experimental parameters that might need to be changed. They are read at runtime and don't require recompilation to change.

## Code parameters
Hard coded at the top of `swapint.cpp` and `swapmotion.cpp`, including parameters like the default configuration file name and the output file name. These shouldn't need to be modified, but if they do, simply change them and recompile.
