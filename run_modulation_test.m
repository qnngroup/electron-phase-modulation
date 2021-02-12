clear all;

% Constants:
c = 299.79245; %speed of light in nm/fs

%Time axis and wavepacket pulse envelope for electron:
t = linspace(-100, 100, 2000); %Time axis in fs
W0 = 1e3; %Central energy in eV
t_electron = 10; %FWHM of electron wavepacket in fs

%Modulator settings:
lambda = 1000; %Wavelength in nm
V_mag = 5; %Modulator potential strength in eV
t_V = 200; %Duration of modulation function (FWHM in fs)

%Demodulator Settings
V_mag_dm = 5;
phase_dm = linspace(0, 8*pi, 100);

%-- Calculations -- 

%Calculate the central frequency of the modulator and its function
omega = 2*pi*c/lambda;
[V_env, garbage] = gaussian_pulse(t, t_V, omega, 0); %Potential function
V = V_mag*V_env; %Construct the actual potential profile

%Calculate the wavepacket envelope
[A, garbage] = gaussian_pulse(t, t_electron, 0, 0);

%Calculate and plot the energy spectrum
[k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V);

figure(1);
plot(W, P_W);
set(gca, 'fontsize', 14);
xlabel('Electron Energy (eV)', 'fontsize', 14);
ylabel('Spectral Amplitude (a.u.)', 'fontsize', 14);

%Now propagate to t_prime and output the real-space wavefunction:
t_prime = linspace(0, 1000, 100);

%Propagate...
[x_center, x, u_out] = propagate_fixed_time(t_prime, W0, k, a_k);

%Plot output wavefunction u_out:
figure(2);
imagesc(x, x_center/1000, abs(u_out).^2);
set(gca, 'fontsize', 14)
xlim([-500, 500]);
shading interp;
colorbar('fontsize', 14);
xlabel('x (nm)', 'fontsize', 14);
ylabel('Propagation Distance (\mum)', 'fontsize', 14);

% -- What about de-modulation? 
physical_constants_normalized;

k0 = sqrt(2*W0);
x_center = 5000;
%t_prime = 1000;
%x_center = t_prime*k0*x0/t0;

%t = linspace(-100, 100, 10000);


[t_center, u_out] = propagate_fixed_space(x_center, t, W0, k, a_k);

for a = 1:length(phase_dm)
  
  [V_env_dm, garbage] = gaussian_pulse(t, t_V, omega, phase_dm(a)); %Potential function
  V_dm = V_mag_dm*V_env_dm; %Construct the actual potential profile
  
  [k_dm, a_k_dm, W_dm, P_W_dm(a, :)] = calc_energy_spec(t, u_out, W0, V_dm);

end

figure(3);
imagesc(phase_dm, W_dm, P_W_dm.');
ylim([950, 1050]); 


