%
% function Mw = m02mw(imag,M0)
% CARL TAPE, 31-Oct-2007
% printed xxx
%
% This function inputs a vector of moments (M0) in units of N-m and outputs
% a column vector of moment magnitudes (Mw).
%
% See Latex notes cmt.pdf.
%
% EXAMPLE (Denali earthquake):
%    M0 = 7.48*1e20; Mw = m02mw(1,M0)
%
% calls xxx
% called by richter.m
%

function Mw = m02mw(imag,M0)

% convert moment tensor from N-m to dyne-cm, since these formulas are
% designed for dyne-cm
M0 = 1e7 * M0(:);

if imag==1
    % Kanamori 1978
    %k = -10.7;
    
    % Kanamori 1977
    %k = -(2/3)*(11.8 - log10(1/2e4));   % for M0 in units of dyne-cm
    %k = -(2/3)*(11.8 - log10(5e-5));    % dyne-cm
    %k = -(2/3)*(16.8 - log10(5));       % dyne-cm
    %k = -11.2 + (2/3)*log10(5);          % dyne-cm
    
    % Kanamori 1977 or 1978
    %Mw = (2/3) * log10(M0) + k;
    
    % Kanamori 1977 (modified form, but exact)
    A = 2/(3*log(10));
    K = 0.2*10^16.8;
    Mw = A*log(M0/K);
    
elseif imag == 2
    % Harvard CMT
    Mw = (2/3) * (log10(M0) - 16.1);
    
else
    error('imag must be 1 or 2.');
end

%=====================================================
