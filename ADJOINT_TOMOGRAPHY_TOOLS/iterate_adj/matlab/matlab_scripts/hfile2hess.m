%
% display_meas.m
% CARL TAPE, 18-April-2008
% printed xxx
%
% /cig/seismo/3D/mt_measure_adj_work/scripts_tomo/matlab/
%
% calls xxx
% called by xxx
%

function H = hfile2hess(hfile,nsrc)

if ~exist(hfile,'file'), error([hfile ' does not exist']); end
[kall,iall,jall,Hij] = textread(hfile,'%f%f%f%f');

n = length(kall);
nhess = (-1 + sqrt(1+8*n))/2;   % quadratic formula
if nhess ~= nsrc, error('number of events must be consistent'); end

% construct the Hessian matrix
H = zeros(nsrc,nsrc);
for i = 1:nsrc
    for j = i:nsrc
        ih = find( and( i == iall, j == jall) );
        H(i,j) = Hij(ih);
    end
end
H = H + H' - diag(diag(H));     % fill the lower triangle

%======================================================
