clear all;

c = 299.79245 %speed of light in nm/fs

t = linspace(-2500, 2500, 100000);
[A, garbage] = gaussian_pulse(t, 200, 0, 0);

lambda = 1000;
omega = 2*pi*c/lambda;

V = sin(omega*t);
%[V, garbage] = gaussian_pulse(t, 200, omega, 0);
V = 10*V;

W0 = 15e3;

[k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V);

figure(1);
plot(W, P_W);

% -- Need to build these up as a function of time and space!
% TODO: Get this commented out and into its own function!
physical_constants_normalized;


k0 = sqrt(2*W0);

k = fftshift(k);
a_k = fftshift(a_k);

k_lin = linspace(k(1), k(end), length(k));
a_k_lin_mag = interp1(k, abs(a_k), k_lin);
a_k_lin_phase = interp1(k, unwrap(angle(a_k)), k_lin);
a_k_lin = a_k_lin_mag.*exp(1j*a_k_lin_phase);

k_lin = fftshift(k_lin);
a_k_lin = fftshift(a_k_lin);

%Now get x in real space
x = linspace(0, 2*pi*x0/(k(2) - k(1)), length(k_lin));

t_out = linspace(0, 4000, 100)/t0;

for a = 1:length(t_out)

  %The wavepacket is moving so we want to keep it in the frame of the simulation...
  %A linear phase ramp gives a space-shift -- we use the central velocity k0 (in normalized units)
  space_shift = exp(-1j*k_lin*k0*t_out(a));
  %space_shift = 1;


  %Phase term that represents propagation in time:
  W_prop = (k_lin.^2/2);
  time_prop = exp(1j*W_prop.*t_out(a));
  u_out(a, :) = fftshift(ifft(a_k_lin.*time_prop.*space_shift));

end



figure(2);
%plot(x, abs(u_out).^2)
imagesc(x, t_out*t0, abs(u_out).^2);
xlim([x(end)/2 - 500, x(end)/2 + 500]);
shading interp;
colorbar();
