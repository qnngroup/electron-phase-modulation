function [E, A] = gaussian_pulse(t, fwhm, wc, phi_ce)
% [E, A] = gaussian_pulse(t, fwhm, wc, phi_ce)
%
% Function for generating a Gaussian pulse form with given intensity FWHM.
% Outputs electric field profile and intensity envelope, both normalized to
% peak of 1.  Units need to be self-consistent between t, fwhm, and wc.
%
% Inputs:
% -----------
%   t --> time 
%   fwhm --> full width at half max of the intensity envelope
%   wc --> central frequency (rad/unit time)
%   phi_ce --> carrier envelope phase offset (rad)
%
% Outputs:
% -------------
%  E --> Electric field profile (peak of 1)
%  A --> Intensity envelope (peak of 1)

%Now determine constant for gaussian envelope from fwhm:
tau = fwhm/(2*sqrt(log(2)));

%Finally, define the pulse:
E_envelope = exp(-t.^2/(2*tau^2));

E = E_envelope.*cos(wc*t + phi_ce);

A = E_envelope.^2;

end

