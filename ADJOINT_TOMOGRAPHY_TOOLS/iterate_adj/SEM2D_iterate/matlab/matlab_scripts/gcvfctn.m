%
% function G = gcvfctn(h, s2, fc, rss0, dof0)
% Carl Tape (Tapio Schneider, ACM 118)
% 09-Nov-2004
%
% GCVFCTN    Evaluate generalized cross-validation function.
%
%    gcvfctn(h, s2, fc, rss0, dof0) is the value of the GCV function
%    for ridge regression with regularization parameter h.
%
%    INPUT:
%       h       regularization parameter
%       s2      squared singular value of the design matrix
%       fc      coefficients fc = U(:, 1:q)'*g
%       rss0    the residual sum of squares of the (non-regularized)
%                   least squares estimate
%       dofo0   the number of residual degrees of freedom of
%                   the (non-regularized) least squares estimate
%
%       U       matrix of left singular vectors
%       g       vector of response variables
%       q       the smaller of the numbers of rows and columns of
%                   the design matrix, q = min(n,p)
%
%    Auxiliary function for GCV.
%
%    Adapted from GCVFUN in Per-Christian Hansen's Regularization Toolbox. 
%
% See Schneider (2001) for details.
%
% calls xxx
% called by gcv.m
%

function G = gcvfctn(h, s2, fc, rss0, dof0)

% SIMILAR TO filter factor for Tikhonov regularization:
% 1 - f = F, where F are the filter factors
f = h^2 ./ (s2 + h^2);

RSS = norm(f.*fc)^2 + rss0;
T2  = (dof0 + sum(f))^2;

G = RSS / T2;

%=============================================================
  