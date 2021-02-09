clear all;

% Constants:
c = 299.79245; %speed of light in nm/fs

%Time axis and wavepacket pulse envelope for electron:
t = linspace(-1000, 1000, 100000); %Time axis in fs
W0 = 1e3; %Central energy in eV
t_electron = 200; %FWHM of electron wavepacket in fs

%Modulator settings:
lambda = 1000; %Wavelength in nm
V_mag = 1; %Modulator potential strength in eV
t_V = 200; %Duration of modulation function (FWHM in fs)


%-- Calculations -- 

%Calculate the central frequency of the modulator and its function
omega = 2*pi*c/lambda;
[V_env, garbage] = gaussian_pulse(t, t_V, omega, 0); %Potential function
V = V_mag*V_env; %Construct the actual potential profile

%Calculate the wavepacket envelope
[A, garbage] = gaussian_pulse(t, 200, 0, 0);

%Calculate and plot the energy spectrum
[k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V);

figure(1);
plot(W, P_W);
set(gca, 'fontsize', 14);
xlabel('Electron Energy (eV)', 'fontsize', 14);
ylabel('Spectral Amplitude (a.u.)', 'fontsize', 14);

##%Now propagate to t_prime and output the real-space wavefunction:
##
##t_prime = linspace(0, 1000, 100);
##
##%Propagate...
##[x_center, x, u_out] = propagate(t_prime, W0, k, a_k);
##
##%Plot output wavefunction u_out:
##figure(2);
##imagesc(x, x_center/1000, abs(u_out).^2);
##set(gca, 'fontsize', 14)
##xlim([-500, 500]);
##shading interp;
##colorbar('fontsize', 14);
##xlabel('x (nm)', 'fontsize', 14);
##ylabel('Propagation Distance (\mum)', 'fontsize', 14);


%PROBLEM: Somehow the time-axis below is stretched incorrectly... I don't know 
%to get it stretched back correctly.  I need to think through how it has been
%converted to this point...

% -- What about de-modulation? 
physical_constants_normalized;
k0 = sqrt(2*W0);
t_prime = 1000;
[x_center, x, u_out] = propagate(t_prime, W0, k, a_k);

%Get the new time-domain envelope
t = -1*t0*fliplr(x/x0/k0);
A = fliplr(u_out);

phase_dm = linspace(0, 4*pi, 50)

 
%Calculate and plot the energy spectrum
for a = 1:length(phase_dm)
  
  V_mag_dm = 1;
  [V_env_dm, garbage] = gaussian_pulse(t, t_V, omega, phase_dm(a)); %Potential function
  V_dm = V_mag_dm*V_env_dm; %Construct the actual potential profile

  [k_dm, a_k_dm, W_dm, P_W_dm(a, :)] = calc_energy_spec(t, A, W0, V_dm);
  
endfor


figure(3);
imagesc(W_dm, phase_dm, P_W_dm);
xlim([990, 1010]); 