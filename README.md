# Electron Phase Modulator Simulations

These functions and scripts are for simulating optically-driven electron phase modulation.  You can find a description the theory that this code is based on along with demonstration usage in [NOTES-AND-DEMOS/electron-phase-modulation-technical-notes-and-demos.ipynb](NOTES-AND-DEMOS/electron-phase-modulation-technical-notes-and-demos.ipynb).  

# Core Functions

## `function [k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V)`

Calculates the energy spectra based on:
 1. A wavefunction at the entrance of a phase plate; and
 2. The oscillating potential energy step due to an oscillating electric field across the phase plate.

It assumes the wavefunction is defined by a pulsed electron wavepacket
with a time-domain complex envelope given by A, and central energy E0.
The slowly-varying envelope approximation is used throughout for these 
calculations. 

For convenience, it outputs both an energy spectrum in eV, and the momentum
spectrum where the k-vector is in the normalized eV units (see below). 

For the calculations below, we use "eV" normalized units for convenience:
 * Energy in eV
 * $\hbar = m = q = 1$

### Inputs 
 * `t` --> time axis (femtoseconds) 
 * `A` --> complex amplitude of wavefunction
 * `W0` --> central energy (in eV)
 * `V` --> time-dependent potential (eV)

### Outputs
 * `k` --> momentum axis (normalized units)
 * `a_k` --> momentum amplitude
 * `W` --> energy axis in eV
 * `P_W` --> energy spectrum
 
## `function [t_center, u_out] = propagate_fixed_space(x_center, t, W0, k, a_k)`

Propagates the wavefunction defined by `k` and `a_k` in time until it is centered around `x_center`.  It then reconstructs the wavefunction in time at that position. 

The output in time is centered such that `t=t_center`.  This is to avoid aliasing and wrapping problems in the fft due to the limited number of data points provided to define the wavefunction.  

    
For the calculations below, we use "eV" normalized units for convenience:
 * Energy in eV
 * $hbar = m = q = 1$

### Inputs
 * `x_center` --> spatial location to evaluate wavefunction (nanometers)
 * `t` --> vector of times giving time range to evaluate wavefunction over
 * `W0` --> central energy in eV
 * `k` --> k-vector (normalized units -- takes output from calc_energy_spec()).
 * `a_k` --> momentum spectrum amplitude (takes output from calc_energy_spec()).

### Outputs
 * `t_center` --> central time of the time-axis in femtoseconds
 * `u_out` --> wavefunction amplitude in time centered around t_center

## `function [x_center, x, u_out] = propagate_fixed_time(t_prime, W0, k, a_k)`

Propagates the wavefunction defined by 'k' and 'a_k' with central energy `W0` forward to time `t_prime`.  It then reconstructs the wavefunction `u_out` in real-space `x` around the position `x_center`.   

Based on the central energy `W0`, it reconstructs the wavefunction in the moving frame.  The output `x` is centered such that `x=0` is at `x_center`.  This is to avoid aliasing and wrapping problems in the fft due to the limited number of data points provided to define the wavefunction. 

For the calculations below, we use "eV" normalized units for convenience:
 * Energy in eV
 * $hbar = m = q = 1$

### Inputs
 * `t_prime` --> time to evaluate the wavefunction (femtoseconds; can be vector)
 * `W0` --> central energy in eV
 * `k` --> k-vector (normalized units -- takes output from `calc_energy_spec()`).
 * `a_k` --> momentum spectrum amplitude (takes output from `calc_energy_spec()`).

### Outputs
 * `x_center` --> central position of the `x`-axis (`x=0`) in nm
 * `x` --> x-axis centered such that `x=0` is at `x_center` in the moving frame of the wavepacket
 * `u_out` --> wavefunction amplitude in real space over `x` in the moving frame for each `t_prime`

# Helper Functions

## `function [E, A] = gaussian_pulse(t, fwhm, wc, phi_ce)`

Function for generating a Gaussian pulse form with given intensity FWHM. Outputs electric field profile and intensity envelope, both normalized to peak of 1.  Units need to be self-consistent between `t`, `fwhm`, and `wc`.  While this was intended to be used to create electric field waveforms, it can generally be used to create Gaussian pulses.  

### Inputs

 * `t` --> time 
 * `fwhm` --> full width at half max of the intensity envelope
 * `wc` --> central frequency (rad/unit time)
 * `phi_ce` --> carrier envelope phase offset (rad)

### Outputs
 * `E` --> Electric field profile (peak of 1)
 * `A` --> Intensity envelope (peak of 1)
 
# Test Script

A test script `run_modulation_test.m` is provided that demonstrates each element of the code for convenience.  