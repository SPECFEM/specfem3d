%
% function [dsub, indmat] = wave2d_dsub(d,covd,nvec)
% Carl Tape, 25-Feb-2010
%
% This function splits a wave2d.f90 model vector into constituent parts.
%
% INPUT
%    d     N x 1 data vector
%    covd  N x 1 diagonal terms of data covariance matrix
%    nvec  S x 1 vector of number of measurments per source
%
% OUTPUT
%    dsub  S x 1 provected data vector
%
% calls xxx
% called by xxx
%

function [dsub, indmat] = wave2d_dsub(d,covd,nvec)

d = d(:);
covd = covd(:);
n = sum(nvec);
nsrc = length(nvec);

if length(d) ~= length(covd), error('d and covd not same lengths'); end
if length(d) ~= n, error('inconsistent number of measurements'); end

% indexing for measurements
cnvec = cumsum(nvec);
indmat = [ [1 ; 1+cnvec(1:end-1)] cnvec ];
disp([indmat nvec]);

dsub = zeros(nsrc,1);
for isrc = 1:nsrc
    inds = [indmat(isrc,1) : indmat(isrc,2)];
    dtemp = d(inds);
    ctemp = covd(inds);
    
    % L2-norm-squared
    %dsub(isrc) = dtemp' * diag(1./ctemp) * dtemp;
    dsub(isrc) = sum( dtemp.^2 ./ ctemp );
end

% check
if any(dsub==0), error(' For at least one source, there is perfect fit.'); end

%=========================================================