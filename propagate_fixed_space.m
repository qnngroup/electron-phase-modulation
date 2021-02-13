function [t_center, u_out] = propagate_fixed_space(x_center, t, W0, k, a_k)
%function [t_center, u_out] = propagate_fixed_space(x_center, t, k, a_k)
%
%Propagates the wavefunction defined by k and a_k with central energy W0 in time until it is centered 
%around x_center.  It then reconstructs the wavefunction in time at that position. 
%
%The output in time is centered such that t=t_center.  This is to 
%avoid aliasing and wrapping problems in the fft due to the limited number of 
%data points provided to define the wavefunction.  
%
% Inputs
% -------
% x_center --> spatial location to evaluate wavefunction (nanometers)
% t --> vector of times giving time range to evaluate wavefunction over
% W0 --> central energy in eV
% k --> k-vector (normalized units -- takes output from calc_energy_spec()).
% a_k --> momentum spectrum amplitude (takes output from calc_energy_spec()).
%
% Outputs
% --------
% t_center --> central time of the time-axis in femtoseconds
% u_out --> wavefunction amplitude in time centered around t_center
%    
%For the calculations below, we use "eV" normalized units for convenience:
% - Energy in eV
% - hbar = m = q = 1

%Load the normalized units constants:
physical_constants_normalized;

%Calculate central momentum:
k0 = sqrt(2*W0);

%Convert to normalized t_prime for calculations below:
t_center = x_center*t0/x0/k0; %central time in femtoseconds
t_prime = (t + t_center)/t0;

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

%Run through each point in t_prime and build up u_out matrix.
%TODO: Can cast this for loop as a matrix operation for speed.
%wait = waitbar(0);
%for a = 1:length(t_prime)

%The wavepacket is moving so we want to keep it in the frame of the simulation...
%A linear phase ramp gives a space-shift -- we use the central velocity k0 (in normalized units).
%Linear phase ramp in k-space is spatial shift:
space_shift = (exp(-1j*k_lin*x_center/x0).')*ones(size(t_prime));

%Phase term that represents propagation in time:
W_prop = (k_lin.^2/2 - W0);
time_prop = exp(1j*(W_prop.')*t_prime);

a_k_lin_mat = (a_k_lin.')*ones(size(t_prime));

u_out_full = ifft(a_k_lin_mat.*time_prop.*space_shift);

u_out = squeeze(u_out_full(1, :));

%waitbar(a/length(t_prime), wait);

%end