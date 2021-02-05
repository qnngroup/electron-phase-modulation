%This is a convenience script for loading the normalized physical
%constants for the TDSE code.  
%
%These units allow you to take energy in eV, and hbar = m = q = 1.
%
%Below we provide some convenient constants -- in particular x0 and t0
%which can be used to convert from nm to the normalized units, and 
%femtoseconds to the normalized units respectively.  

hb = 1.05457173e-34; % hbar
me = 9.10938291e-31; % electron mass
e0 = 1.60217657e-19; % elementary charge
eps0 = 8.85418782e-21; % permittivity of free-space (F/nm)

x0 = hb/sqrt(me*e0)/(1e-9); % unit of length ~0.276 nm
t0 = hb/e0/(1e-15); % unit of time ~0.658 fs (energy in eV)
c = 299.792458 * (t0 / x0); % speed of light in eV units
