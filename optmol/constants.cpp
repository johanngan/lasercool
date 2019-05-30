#include "constants.hpp"

const double HBAR = 1.054571800e-34;
const double SPEED_OF_LIGHT = 299792458;
const double K_BOLTZMANN = 1.380649e-23;
// Relative mass of neutral Be = 9.0121831 u
// u = 1/12 g/mol
// mol = 6.02214076e23
// electron mass = 9.10938e-31 kg
const double MASS_BE_PLUS = 1.2461792e-27;
// decay_rate = 2*kB*T_doppler/hbar
// T_doppler = 0.47 mK for Be+
const double DECAY_RATE_BE_PLUS = 1.23e8;
// 1/wavelength = 31928.7 1/cm
// k = 2*pi/wavelength
const double WAVENUMBER_BE_2S2P = 1.00307e7;

const double MASS_RB = 1.411e-25;
const double DECAY_RATE_RB = 1/27e-9;
// 1/wavelength = 12816.545 1/cm for 2P(3/2)->2S
const double WAVENUMBER_RB = 8.0528727e6;
