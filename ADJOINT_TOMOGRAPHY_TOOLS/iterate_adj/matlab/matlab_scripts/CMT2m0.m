%
% function M0 = CMT2m0(im0,M)
% CARL TAPE, 02-Feb-2007
% printed xxx
%
% This function converts from a (CMT) moment tensor to scalar seismic
% moment.  The units of M0 are the same as the elements of Mij, which
% should be Newton-meter (N-m), although it does not affect the function
% here.
%
% See Ekstrom email (11-Oct-2006) and corresponding Latex notes.
%
% moment tensor M = [Mrr Mtt Mpp Mrt Mrp Mtp]
%
% calls xxx
% called by richter.m, denaliquake.m
%

function M0 = CMT2m0(im0,M)

% make sure M is N by 6
[a,b] = size(M); if b ~= 6, M = M'; end
[a,b] = size(M); if b ~= 6, error('dimension of M must be N by 6'); end

ncmt = a;
Mrr = M(:,1); Mtt = M(:,2); Mpp = M(:,3);
Mrt = M(:,4); Mrp = M(:,5); Mtp = M(:,6);

M0 = zeros(ncmt,1);
for ii=1:ncmt
    
    % if you need to compute the eigenvalues
    if or(im0==2,im0==3)
        % convention: r (up), theta (south), phi (east)
        Mcmt = zeros(3,3);
        Mcmt = [ Mrr(ii) Mrt(ii) Mrp(ii) ;
                 Mrt(ii) Mtt(ii) Mtp(ii) ;
                 Mrp(ii) Mtp(ii) Mpp(ii) ];

        [V, D] = eig(Mcmt);
        lams = diag(D)';
        isign = sign(lams);
        %(inv(V)*Mcmt*V - D) / norm(Mcmt)

        % adjust the order of eigenvectors
        % [lamsort, isort] = sort(abs(lams),'descend');
        % Vsort = V(:,isort)
        % Dsort = diag(lamsort.*isign(isort))
        % (inv(Vsort)*Mcmt*Vsort - Dsort) / norm(Mcmt)  % check
    end

    % formula to convert M to M0
    % (1) Dahlen and Tromp, eq. 5.91
    % (2) Harvard CMT way -- double-couple moment
    % (3) Caltech way
    switch im0
        case 1, M0(ii) = 1/sqrt(2) * sqrt( Mrr(ii)^2 + Mtt(ii)^2 + Mpp(ii)^2 ...
                + 2*(Mrt(ii)*Mrt(ii) + Mrp(ii)*Mrp(ii) + Mtp(ii)*Mtp(ii) ) );
        case 2, M0(ii) = (abs(min(lams)) + abs(max(lams)) ) / 2;
        case 3, M0(ii) = max(abs(lams));
    end
end

%=====================================================
