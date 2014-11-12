%
% function rss_vec = ocv_carl(g, A, lamvec)
% Carl Tape, 30-March-2006
%
% Copied from gcv_carl.m on 27-March-2006.
%
% Returns the ordinary cross-validation (OCV) function corresponding to a
% set of input regularization (or damping) parameters.
%
% TWO ALGORITHMS ARE SHOWN:
%   (1) Brute force, which involves ndata inversions per lambda
%   (2) Elegant formalisum, which involves one inversion per lambda
%       --> See Latex notes
%       /home/carltape/classes/acm118/2006_handouts/hw3/hw3sol_2006_prob3.pdf
%
% Using some GPS data, I checked that these approaches are identical.
%
% calls xxx
% called by ridge_carl.m
% 

function rss_vec = ocv_carl(d, A, lamvec)

% Size of inputs
[ndata, nparm]  = size(A); 
numlam          = length(lamvec);

if (min(lamvec) < 0)
    error('Impossible regularization parameter lambda.')
end

rss_vec = zeros(numlam,1);

% loop over regularization parameters
for ii=1:numlam
    lam = lamvec(ii); 
    disp([' ii = ' num2str(ii) '/' num2str(numlam) ', lam = ' num2str(lam)]);
    
    if 1==1
        H = A*inv(A'*A + lam^2*eye(nparm))*A';
        dhat = H*d;
        
        % OCV residual
        res = (d - dhat) ./ (1 - diag(H));
        
        % sum the residuals
        rss_vec(ii) = sum(res.^2);
        
    else
        % loop over datapoints
        for jj=1:ndata
            %disp([' jj = ' num2str(jj) '/' num2str(ndata) ]);

            % indices for which you compute the model parameters
            switch jj
                case 1,     einds = [2:ndata];
                case ndata, einds = [1:ndata-1];
                otherwise,  einds = [1:jj-1 jj+1:ndata];
            end

            % indices to estimate the RSS
            oinds = jj;

            % reduced matrices
            X = A(einds,:);
            g = d(einds);

            % note: regularization matrix is identity matrix
            f_h = inv(X'*X + lam^2*eye(nparm))*X'*g;

            % estimate the model at the datapoints NOT used
            res = A(oinds,:) * f_h - d(oinds);

            % sum the residuals
            rss_vec(ii) = rss_vec(ii) + sum(res.^2);
        end
    end
    
end

rss_vec = rss_vec / ndata;  % normalization

figure; loglog(lamvec, rss_vec, '.'); grid on;
xlabel(' Regularization parameter, log (\lambda)');
ylabel(' RSS at datapoints NOT used in computing model');
   
%======================================================
