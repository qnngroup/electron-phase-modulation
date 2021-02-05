clear all;

c = 299.79245 %speed of light in nm/fs

t = linspace(-500, 500, 50000);
[A, garbage] = gaussian_pulse(t, 250, 0, 0);

lambda = 1000;
omega = 2*pi*c/lambda;

V = sin(omega*t);
%[V, garbage] = gaussian_pulse(t, 200, omega, 0);
V = V*1;

W0 = 1e3;

[k, a_k, W, P_W] = calc_energy_spec(t, A, W0, V);

figure(1);
plot(W, P_W);

% -- Need to build these up as a function of time and space!
% TODO: Get this commented out and into its own function!
physical_constants_normalized;
t_out = 2750/t0;
k0 = sqrt(2*W0);
space_shift = exp(+1j*k*k0*t_out);
W_prop = (k.^2/2);
time_prop = exp(-1j*W_prop.*t_out);
u_out = ifft(a_k.*time_prop.*space_shift);

%Now get x in real space
x = linspace(0, 2*pi*x0/(k(1) - k(end)), length(k));

figure(2);
plot(x, abs(u_out))
