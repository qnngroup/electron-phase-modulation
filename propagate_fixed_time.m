function [x_center, x, u_out] = propagate_fixed_time(t_prime, W0, k, a_k)
%function [x, u_out] = propagate(x_center, k, a_k)
%
%Propagates the wavefunction defined by k and a_k forward to time t_prime.  It
%then reconstructs the wavefunction in real-space.  
%
%Based on the central energy W0, it reconstructs the wavefunction in the moving
%frame.  The output x is centered such that x=0 is at x_center.  This is to 
%avoid aliasing and wrapping problems in the fft due to the limited number of 
%data points provided to define the wavefunction.  
%
% Inputs
% -------
% t_prime --> time to evaluate the wavefunction (femtoseconds; can be vector)
% W0 --> central energy in eV
% k --> k-vector (normalized units -- takes output from calc_energy_spec()).
% a_k --> momentum spectrum amplitude (takes output from calc_energy_spec()).
%
% Outputs
% --------
% x_center --> central position of the x-axis (x=0) in nm
% x --> x-axis centered such that x=0 is at x_center in the moving frame of the wavepacket
% u_out --> wavefunction amplitude in real space over x in the moving frame for each t_prime
%    
%For the calculations below, we use "eV" normalized units for convenience:
% - Energy in eV
% - hbar = m = q = 1

%Load the normalized units constants:
physical_constants_normalized;

%Calculate central momentum:
k0 = sqrt(2*W0);

%Convert to normalized t_prime for calculations below:
t_prime = t_prime/t0;

%Prep the input data for interpolation:
k = fftshift(k);
a_k = fftshift(a_k);

%This needs to be linearized -- we calculated in time/energy at first which gives
%a quadratic momentum axis -- for FT to work using fft we need a linearized k-axis.
k_lin = linspace(k(1), k(end), length(k));
a_k_lin_mag = interp1(k, abs(a_k), k_lin);
a_k_lin_phase = interp1(k, unwrap(angle(a_k)), k_lin);
a_k_lin = a_k_lin_mag.*exp(1j*a_k_lin_phase);

%Shift it back...
k_lin = fftshift(k_lin);
a_k_lin = fftshift(a_k_lin);

%Now get x in real space
x = linspace(0, 2*pi*x0/(k(2) - k(1)), length(k_lin));
x = x - x(end)/2;

%Center position of wavepacket window
x_center = t_prime*k0*x0; %in nm

%Run through each point in t_prime and build up u_out matrix.
%TODO: Can cast this for loop as a matrix operation for speed.
for a = 1:length(t_prime)

  %The wavepacket is moving so we want to keep it in the frame of the simulation...
  %A linear phase ramp gives a space-shift -- we use the central velocity k0 (in normalized units).
  %Linear phase ramp in k-space is spatial shift:
  space_shift = exp(-1j*k_lin*k0*t_prime(a));

  %Phase term that represents propagation in time:
  W_prop = (k_lin.^2/2 - W0);
  time_prop = exp(1j*W_prop.*t_prime(a));
  
  
  u_out(a, :) = fftshift(ifft(a_k_lin.*time_prop.*space_shift));

end
