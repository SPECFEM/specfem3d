% CRITICAL_TIMESTEP	Computes the critical timestep in 3D, assuming:
%                       purely elastic medium
%                       the leapfrog time scheme
%                       a cube element (same size for all faces)
%
% SYNTAX	dtc = critical_timestep(csp,h,ngll)
%
% INPUT		cp	P wave speed (in m/s)
%		h	element size (in m)
%		ngll	number of GLL nodes per spectral element edge (2 to 20)
%
% OUTPUT	dtc	critical timestep
%
function dtc = critical_timestep(csp,h,ngll)

DIM = 3; % dimension

if ngll<2 | ngll>20, error('ngll must be from 2 to 20'); end

% tabulated values of critical frequency (non-dimensional, 1D)
% Omega_max(p) = Omega_max(ngll-1)
Omega_max = [2.0000000e+00 4.8989795e+00 8.6203822e+00 1.3540623e+01 1.9797952e+01 2.7378050e+01 ...
     3.6256848e+01 4.6421894e+01 5.7867306e+01 7.0590158e+01 8.4588883e+01 9.9862585e+01 ...
     1.1641072e+02 1.3423295e+02 1.5332903e+02 1.7369883e+02 1.9534221e+02 2.1825912e+02 2.4244948e+02];

% stability factor por leapfrog timescheme
C = 2;

% critical time step, 
% assumes a cube element
dtc = C*h/csp/sqrt(DIM)/Omega_max(ngll-1);
