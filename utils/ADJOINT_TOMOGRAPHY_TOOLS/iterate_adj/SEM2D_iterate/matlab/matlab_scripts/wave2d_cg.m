%
% function [m,pk] = wave2d_cg(m0,gk,C,chi_k_val,mu_val,istep,imaketest,gfile,chitfile)
% Carl Tape, 21-Jan-2010
%
% This emulates the CG algorithm in wave2d.f90, but assumes a full
% covariance matrix instead of a diagonal covariance matrix.
%
% calls xxx
% called by xxx
%

function [m,pk] = wave2d_cg(m0,gk,C,chi_k_val,mu_val,istep,imaketest,gfile,chitfile)

disp('------CG ALGORITHM---------');

whos m0 gk C

nmod = length(m0);
mt = zeros(nmod,1);
mk = zeros(nmod,1);

if istep <= 1
    beta_val = 0.0;
    p0 = zeros(nmod,1);
else
    % load p0 and g0, the PREVIOUS gradient vectors (pk and gk)
    if ~exist(gfile,'file')
        disp(gfile);
        error('gfile does not exist');
    end
    [g0,p0] = textread(gfile,'%f%f');
    
    beta_val = ( (gk - g0)' * C * gk ) / ( g0' * C * g0 );
    if isinf(beta_val), error('beta_val is infinity'); end
end
pk = -C * gk + beta_val * p0;
lam_t_val_bot = sum( gk .* pk );        % gk is hat, pk is non-hat
lam_t_val = 2.0*(mu_val - chi_k_val) / lam_t_val_bot;

% check
if 0==1
    disp('checking values in the CG algorithm');
    vals = [mu_val lam_t_val_bot lam_t_val beta_val chi_k_val];
    %write(19,'(5e16.8)') mu_val, lam_t_val_bot, lam_t_val, beta_val, chi_k_val
    dirbase = '/home/carltape/ADJOINT_TOMO/iterate_adj/';
    vals0 = load([dirbase 'SEM2D_iterate_OUTPUT/run_9150/cg_test_vals.dat']);
    for ii=1:length(vals), vals0(ii), vals(ii), end
    error('testing here');
end

if imaketest==1
    % test model
    mt = m0 + lam_t_val*pk;
    
else
    % load chi for test model
    if ~exist(chitfile,'file'), error('chitfile does not exist'); end
    chi_t_val = load(chitfile);
    
    % a quadratic fit requires at least 5 values
    xx1 = 0.0;
    xx2 = lam_t_val;
    yy1 = chi_k_val;
    yy2 = chi_t_val;
    g1  = sum(gk .* pk);
    %g1  = sum(g0 .* pk);
    
    % coefficients of the quadratic polynomial (ax^2 + bx + c)
    Pa = ((yy2 - yy1) - g1*(xx2 - xx1)) / (xx2^2 - xx1^2);
    Pb = g1;
    %Pc = yy1 - Pa*xx1^2 - Pb*xx1;
    
    % get the analytical minimum (the vertex)
    if (Pa ~= 0.0)
        xmin = -Pb / (2.0*Pa);
    else
        error('check the quadratic input polynomial');
    end
    
    % compute updated model
    lam_0_val = xmin;
    mk = m0 + lam_0_val * pk;
end

% assign output model
if imaketest == 1
    m = mt;
else
    m = mk;
end

%=============================================================
