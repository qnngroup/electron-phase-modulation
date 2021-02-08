function [k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V)
%function [k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V)
%
%Calculates the energy spectra based on:
% (1) A wavefunction at the entrance of a phase plate; and
% (2) The oscillating potential energy step due to an oscillating electric field
%     across the phase plate.
%
%It assumes the wavefunction is defined by a pulsed electron wavepacket
%with a time-domain complex envelope given by A, and central energy E0.
%The slowly-varying envelope approximation is used throughout for these 
%calculations. 
%
%For convenience, it outputs both an energy spectrum in eV, and the momentum
%spectrum where the k-vector is in the normalized eV units (see below).  
%
% Inputs
% -------
% t --> time axis (femtoseconds) 
% A --> complex amplitude of wavefunction
% W0 --> central energy (in eV)
% V --> time-dependent potential (eV)
%
% Outputs
% --------
% k --> momentum axis (normalized units)
% a_k --> momentum amplitude
% W --> energy axis in eV
% P_W --> energy spectrum
    
%For the calculations below, we use "eV" normalized units for convenience:
% - Energy in eV
% - hbar = m = q = 1
physical_constants_normalized

%Some initial setup...
t = t/t0; %convert time to normalized units
k0 = sqrt(2*W0); %Central momentum in normalized units
dt = t(2) - t(1);
Nt = length(t);
								
%Calculate the phase modulation phi_mod due to V_t
phi_mod = cumtrapz(t, V);

%Setup energy/momentum values
T = t(end) - t(1);
W = (2*pi/T)*fftshift((-1*floor(Nt/2) - mod(Nt, 2)):(floor(Nt/2)-1)) + W0;
k = sqrt(2*W);

%Calculate the spectral amplitude using ifft:
a_k = ((k + k0)/(2*sqrt(2*pi))).*ifft(fftshift(A.*exp(-1j*phi_mod)));

%For the Jacobian
dk = diff(k);
dk(end + 1) = dk(end);
dW = W(2) - W(1);

%Finally the energy spectrum
P_W = abs(a_k).^2.*dk./dW;

W = ifftshift(W);
P_W = ifftshift(P_W);
