%
% function 
% Carl Tape, 25-Jan-2010
%
% Compute the Hessian for the source subspace projection method.
%
% INPUT
%    Gstr       gradient w.r.t. structure parameters
%    Gsrc       gradient w.r.t. source parameters 
%    Cstr       covariance matrix or diagonal part (as a vector)
%    csrc       diagonal covariance matrix (as a vector)
%
% calls xxx
% called by xxx
%

function [Hstr,Hsrc] = wave2d_compute_hessian(Gstr,Gsrc,Cstr,csrc)

[nsrc,nstr] = size(Gstr);
[a,b] = size(Cstr);
if or(a==1,b==1)
    icdiag = 1;
    Cstr = Cstr(:)';      % ensure that C is a row vector
else
    icdiag = 0;
end

% Hessian for structure parameters (from event kernels)
Hstr = zeros(nsrc,nsrc);
if icdiag == 1
    disp('constructing the Hessian with diagonal Cm for structural part');
    for ii = 1:nsrc
        for jj = 1:nsrc
            Hstr(ii,jj) = dot(Gstr(ii,:), Cstr.*Gstr(jj,:));
        end
    end
    
else
    disp('constructing the Hessian with full Cm for structural part');
    Hstr = Gstr * Cstr * Gstr';
end

% Hessian for source parameters
Hsrc = zeros(nsrc,nsrc);
Hsrc = Gsrc * diag(csrc) * Gsrc';   % csrc is Msrc x 1

disp('norm(Hstr), norm(Hsrc):');
norm(Hstr), norm(Hsrc)

%=========================================================

if 0==1
    S = 5;              % number of sources
    Mstr = 20;
    Msrc = S*3;
    inds_str = [1:Mstr]';
    inds_src = [1:Msrc]';
    
    % fill G
    Gstr = rand(S,Mstr);
    Gsrc = zeros(S,Msrc);
    %Gsrc = repmat(diag(rand(S,1)),1,Msrc1+Msrc2);
    Cstr  = rand(Mstr,Mstr);        % full matrix
    csrc  = rand(Msrc,1);           % vector representing diagonal matrix
    
    % testing
    Hstr = wave2d_compute_hessian(Gstr,Gsrc,Cstr,csrc)
    Hstr = wave2d_compute_hessian(Gstr,Gsrc,diag(diag(Cstr)),csrc)
    Hstr = wave2d_compute_hessian(Gstr,Gsrc,diag(Cstr),csrc)
end

%=========================================================
