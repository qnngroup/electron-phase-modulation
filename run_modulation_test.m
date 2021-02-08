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


%Now propagate to t_prime and output the real-space wavefunction:
t_prime = linspace(0, 1000, 100);

%Propagate...
[x_center, x, u_out] = propagate(t_prime, W0, k, a_k);

%Plot output wavefunction u_out:
figure(2);
imagesc(x, x_center/1000, abs(u_out).^2);
set(gca, 'fontsize', 14)
xlim([-500, 500]);
shading interp;
colorbar('fontsize', 14);
xlabel('x (nm)', 'fontsize', 14);
ylabel('Propagation Distance (\mum)', 'fontsize', 14);
