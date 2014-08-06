%
% wave2d_cg_run.m
% Carl Tape, 21-Jan-2010
%
% This program computes an updated model using a conjugate gradient
% algorithm.  This is essentially a TEST PROGRAM in anticipation of
% performing the inversion in Matlab OUTSIDE of wave2d.f90.
%
% calls xxx
% called by xxx
%

format long
format compact
close all
clear

% add path to additional matlab scripts
path([pwd '/matlab_scripts'],path);

colors;
stfm = '%4.4i';
%stpars = {'B = ln(beta/beta0)','Ts = ts = ts0','Xs = xs - xs0','Ys = ys - ys0'};
stpars = {'B = ln(beta/beta0)','ts','xs','ys'};
npts = 100;

%---------------------------------------------------------
% USER INPUT

% base directory
%dirbase = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/';
dirbase = '/home/carltape/ADJOINT_TOMO/iterate_adj/';
if ~exist(dirbase,'dir'), error(['dirbase ' dirbase ' does not exist']); end

% directory with Matlab parameterization for full covariance matrix
im = 1;     % KEY: index for model (im=1 for 50km, im=2 for 25km)
stim = sprintf('%3.3i',im);
dirrand = [dirbase '/SEM2D_iterate_INPUT/random_fields/model_' stim '/'];
if ~exist(dirrand,'dir'), error(['dirrand ' dirrand ' does not exist']); end

% KEY
%parms = [4000 31];
%parms = [6000 31];
parms = [6000 8];

% icg = 1: emulate wave2d.f90 algorithm (diagonal Cm)
% icg = 2: emulate wave2d.f90 algorithm, but use 32 x 32 cells (diagonal Cm)
% icg = 3: CG algorithm with full Cm and 32 x 32 cells
% icg = 4: source subspace algorithm with full Cm and 32 x 32 cells
icg = 4;
iwrite = 1;
iplotmod = 1;
iplotker = 0;
sticg = sprintf('%2.2i',icg);

nsrc = 25;

%---------------------------------------------------------

% whether to use smoothed event kernels
if any(icg == [3 4]), ismooth = 0; else ismooth = 1; end

irun0 = parms(1);
istep = parms(2);
imaketest = ~mod(istep,2);

if imaketest==1
    im0 = istep;
    imt = istep + 1;
else
    im0 = istep - 1;
    imt = istep;
end
imk = im0+2;
imp = im0-2;
disp(sprintf(' im0 = %i, imt = %i, imk = %i, imaketest = %i',im0,imt,imk,imaketest));

% base run directory
stirun0 = sprintf(stfm,irun0);
dir0 = [dirbase 'SEM2D_iterate_OUTPUT/run_' stirun0 '/'];
if ~exist(dir0,'dir'), error(['dir0 ' dir0 ' does not exist']); end

% model directories
dirR   = [dir0 'READ_IN/'];
if ~exist(dirR,'dir'), mkdir(dirR); end
stimp = sprintf(stfm,imp); stim0 = sprintf(stfm,im0);
stimt = sprintf(stfm,imt); stimk = sprintf(stfm,imk);
dirmp  = [dirR 'model_m' stimp '/'];
dirm0  = [dirR 'model_m' stim0 '/'];
dirmt  = [dirR 'model_m' stimt '/'];
dirmk  = [dirR 'model_m' stimk '/'];

% for the first step, read from the base run directory
if istep <= 1, dirm0 = dir0; end
if istep <= 3, dirmp = dir0; end
if ~exist(dirm0,'dir'), error(['dirm0 ' dirm0 ' does not exist']); end

% assign the input and output directories
idir1 = dirm0;
idir2 = dirmp;
if imaketest==1
    imo = imt;
    odir = dirmt;
    if icg==4, imo = imk; odir = dirmk; end
else
    imo = imk;
    odir  = dirmk;
    if ~exist(dirmt,'dir'), error(['dirmt ' dirmt ' does not exist']); end
end
stimo = sprintf(stfm,imo);
if ~exist(idir1,'dir'), error(['idir1 ' idir1 ' does not exist']); end
if ~exist(idir2,'dir'), error(['idir2 ' idir2 ' does not exist']); end

if iwrite==1
    if ~exist(odir,'dir'), mkdir(odir); disp(['making ' odir]); end
end

%---------------------------------------------------------
% load files

% read parameter file
[stnames,stvals] = textread([idir1 'wave2d_constants.dat'],'%s%f');
for ii=1:length(stnames), eval([stnames{ii} ' = stvals(' num2str(ii) ')']); end

% constants for model covariance
[stnames,stvals] = textread([idir1 'scaling_values_covm.dat'],'%s%f');
%for ii=1:length(stnames), eval([stnames{ii} ' = stvals(' num2str(ii) ')']); end
covm_weight_parts = stvals(1:4);
covg_weight_parts = stvals(6:9);
ugsq_str = stvals(11);
ugsq_ts = stvals(12);
ugsq_xs = stvals(13);
ugsq_ys = stvals(14);
dnevent_run = stvals(15);
coverage_str = stvals(16);
coverage_src = stvals(17);

% constants for data covariance
[stnames,stvals] = textread([idir1 'scaling_values_covd.dat'],'%s%f');
ievent_min = stvals(1);
ievent_max = stvals(2);
nevent_run = stvals(3);
nrec       = stvals(4);             % number of stations
ncomp      = stvals(5);             % number of components
nmeas_run  = stvals(6);             % number of measurements
sigma_DT   = stvals(7);             % std errors added to DT measurements

nparm_src_inv = sum([INV_SOURCE_T 2*INV_SOURCE_X ]) * nevent_run;

% event indices
einds = [ievent_min:ievent_max];

% %   double precision, parameter :: DENSITY           = 2.60e3 ! kg/m^3
% %   double precision, parameter :: INCOMPRESSIBILITY = 5.20e10 ! Pa
% %   double precision, parameter :: RIGIDITY          = 2.66e10 ! Pa
% rho0 = 2.60e3;
% kap0 = 5.20e10;
% mu0 = 2.66e10;
% beta0 = sqrt( mu0/rho0 )
% bulk0 = sqrt( kap0 /rho0 )

% indexing
m_inds = load([idir1 'm_inds.dat']);
nmod = m_inds(4,2);
nlocal = m_inds(1,2);
indB = m_inds(1,1):m_inds(1,2);
indT = m_inds(2,1):m_inds(2,2);
indX = m_inds(3,1):m_inds(3,2);
indY = m_inds(4,1):m_inds(4,2);
indTXY = m_inds(2,1):m_inds(4,2);
nmod_str = nlocal;
nmod_src = nmod - nmod_str;

% local mesh
% write(15,'(6i8,4e18.8)') k, ispec, i, j, iglob, valence(iglob), &
%    x(iglob), z(iglob), da_local(i,j,ispec), da_global(iglob)
[i1,i2,i3,i4,i5,i6,xg,yg,da_local_vec,d10] = textread([idir1 'local_mesh.dat'],'%f%f%f%f%f%f%f%f%f%f');
xmin = min(xg); xmax = max(xg);
ymin = min(yg); ymax = max(yg);
Atot = (xmax-xmin)*(ymax-ymin);
% check
sum(da_local_vec), (xmax-xmin)*(ymax-ymin)

% to match notation for full covariance matrix
Avec = sqrt( da_local_vec/Atot );
sum(Avec.^2)         % check

% plotting mesh
xp = linspace(xmin,xmax,npts);
yp = linspace(ymin,ymax,npts);
[Xp,Yp] = meshgrid(xp,yp);

% load current model
%[m_str_lon,m_str_lat,m_str_kappa,m_str_mu,m_str_rho,m_str_B] = ...
%    textread([idir1 'structure_syn.dat'],'%f%f%f%f%f%f%f');

% prior, initial, and target models
[mprior, m0_initial, mtarget] = textread([dir0 'prior_initial_target_models.dat'],'%f%f%f');

% CURRENT model vector
%write(19,'(4e16.6)') m0(i), mt(i),  m0_vec(i), mt_vec(i)
[m0, mt0, m0_vec, mt_vec0] = textread([idir1 'cg_model_vectors.dat'],'%f%f%f%f');
%[m0_B, m0_T, m0_X, m0_Y] = wave2d_splitm(m0,m_inds);
%[m0_vec_B, m0_vec_T, m0_vec_X, m0_vec_Y] = wave2d_splitm(m0_vec,m_inds);

% data covariance
% nmeas_run = nevent_run * nrec * NCOMP
%  double precision, parameter :: SIGMA_DT = 0.10
%cov_data(:) = SIGMA_DT * SIGMA_DT * nmeas_run
covd_diag = load([idir1 'cov_data_diagonal.dat']);
covd_const = sigma_DT^2 * nmeas_run;
% check
min(covd_diag), max(covd_diag), covd_const

% data vector
% chi_data(ievent,irec,icomp,1) = (tshift_xc_pert )**2 / cov_data(imeasure)
% write(19,'(3i8,1e20.10)') ievent, irec, icomp, chi_data(ievent,irec,icomp,1)
chi_data_all = load([idir1 'chi_data_all.dat']);
chi_data_vec = chi_data_all(:,4);
chi_data_stop = load([idir1 'chi_data_stop.dat']);

% measurement vector
%measure_vec(imeasure,1) = tshift_xc_pert
%measure_vec(imeasure,2) = tshift_xc
%measure_vec(imeasure,3) = dlnA_pert
%measure_vec(imeasure,4) = dlnA
%measure_vec(imeasure,5) = 0.5 * DT*sum( adj_syn(:,icomp,irec)**2 )
meas_all = load([idir1 'measure_vec.dat']);
DTvec = meas_all(:,1);
data_norm = sum(DTvec.^2 / covd_const);
% check
data_norm0 = load([idir1 'data_norm.dat']);
model_norm0 = load([idir1 'model_norm.dat']);
chi0 = load([idir1 'chi.dat'])
data_norm0, data_norm, sum(chi_data_vec), DTvec' * diag(1./covd_diag) * DTvec

% vector of number of measurement windows per event
Ns_vec = nrec * ones(nsrc,1);   % KEY: same measurements per event
Nsi_vec = zeros(nmeas_run,1);
cNs = cumsum(Ns_vec);
Ns_inds = [ [1 ; cNs(1:end-1)+1] cNs ];

% KEY: compute projected data vector
[dnorm2, indmat] = wave2d_dsub(DTvec,covd_diag,Ns_vec);
whos chi_data_vec dnorm2
data_norm0, data_norm, sum(chi_data_vec), sum(dnorm2)

% compute variance reduction for EVERY SINGLE WINDOW
if and(imaketest==1, istep >= 2)
    % load measurements for PREVIOUS model
    meas_allp = load([idir2 'measure_vec.dat']);
    DTvecp = meas_allp(:,1);
    
    [VR_i, VR_s, VR_tot] = wave2d_VR(DTvecp,DTvec,sigma_DT,covd_diag,Ns_inds);
end

% KEY: stopping criterion
if and(imaketest==1, istep >= 2)
    % previous chi value
    chip = load([dirmp 'chi.dat']);
    disp('----------------');
    chip, chi0, log( chip / chi0 ), VAR_RED_MIN
    if log( chip / chi0 ) < VAR_RED_MIN
        error('VR stopping criteria has been met');
    end
end

%----------------------------------

% sigma values
% write(19,'(2e20.10)') sigma_beta, m_scale_str(1)
% write(19,'(2e20.10)') sigma_ts, m_scale_src(1)
% write(19,'(2e20.10)') sigma_xs, m_scale_src(2)
% write(19,'(2e20.10)') sigma_zs, m_scale_src(3)
sigma_all  = load([idir1 'sigma_values.dat']);
sigma_beta = sigma_all(1,1);
sigma_ts   = sigma_all(2,1);
sigma_xs   = sigma_all(3,1);
sigma_zs   = sigma_all(4,1);

% reference values
[alpha0, beta0, rho0, bulk0, kappa0, mu0] = ...
    textread([idir1 'reference_values.dat'],'%f%f%f%f%f%f');

% load kernels
kernels_all = zeros(nlocal,nevent_run);
for ii=1:nevent_run
    ievent = einds(ii)
    dirK = [idir1 'event_' sprintf('%3.3i',ievent) '/'];
    kernel = load([dirK 'kernel_basis_smooth']);
    if ismooth == 0
        Kbeta = kernel(:,4); cfac = 0.2;
    else
        Kbeta = kernel(:,3); cfac = 0.8;
    end
    kernels_all(:,ii) = Kbeta;
    
    if and(iplotker==1,icg==1)
        Zp = griddata(xg,yg,Kbeta,Xp,Yp,'nearest');
        figure; imagesc(Zp); set(gca,'ydir','normal');
        axis equal, axis tight
        colormap(seis); caxis([-1 1]*max(abs(Kbeta))*cfac);
        title(sprintf('Event kernel %i/%i',ii,nevent_run));
    end
end
summed_ker = sum(kernels_all,2);

% source gradient
source_gradient = load([idir1 'source_gradient.dat']);
%source_gradient = zeros(nmod_src,1);

%========================================================================
% ONE STEP OF A CONJUGATE GRADIENT ALGORITHM

% files for CG algorithm
gfilemp = ['cg' sticg '_grad_m' stimp '.dat'];
gfilem0 = ['cg' sticg '_grad_m' stim0 '.dat'];
chitfile = [dirmt 'chi.dat'];
%---------------------------------------------------------
% icg = 1 ==> emulate CG variables in wave2d.f90
if icg==1
    
    if 1==1
        cov_model = zeros(nmod,1);
        cov_model(indB) = sigma_beta^2 * Atot ./ da_local_vec * coverage_str;
        cov_model(indT) = sigma_ts^2 * dnevent_run * coverage_src;
        cov_model(indX) = sigma_xs^2 * dnevent_run * coverage_src;
        cov_model(indY) = sigma_zs^2 * dnevent_run * coverage_src;
        
        cov_imetric = zeros(nmod,1);
        cov_imetric(indB) = cov_model(indB) / covg_weight_parts(1);
        cov_imetric(indT) = cov_model(indT) / covg_weight_parts(2);
        cov_imetric(indX) = cov_model(indX) / covg_weight_parts(3);
        cov_imetric(indY) = cov_model(indY) / covg_weight_parts(4);
        
    else
        % ONLY FOR TESTING THE COMPARISON
        cov_model = zeros(nmod,1);
        cov_model(indB) = sigma_beta^2;
        cov_model(indT) = sigma_ts^2;
        cov_model(indX) = sigma_xs^2;
        cov_model(indY) = sigma_zs^2;
        
        cov_weight = zeros(nmod,1);
        %cov_weight(indB) = coverage_str / covg_weight_parts(1) * Atot ./ da_local_vec;
        cov_weight(indB) = coverage_str / covg_weight_parts(1) ./ Avec.^2;
        cov_weight(indT) = dnevent_run * coverage_src / covg_weight_parts(2);
        cov_weight(indX) = dnevent_run * coverage_src / covg_weight_parts(3);
        cov_weight(indY) = dnevent_run * coverage_src / covg_weight_parts(4);
        
        cov_imetric = zeros(nmod,1);
        cov_imetric(indB) = cov_model(indB) .* cov_weight(indB);
        cov_imetric(indT) = cov_model(indT) .* cov_weight(indT);
        cov_imetric(indX) = cov_model(indX) .* cov_weight(indX);
        cov_imetric(indY) = cov_model(indY) .* cov_weight(indY);
        
        wave2d_diff_vec(cov_imetric0,cov_imetric,m_inds,stpars,0);
        break
    end
    
    % cov_model -- TO CHECK
    [cov_model0,cov_imetric0] = textread([idir1 'cov_model_diagonal.dat'],'%f%f');
    % check
    disp('comparing cov_model and cov_imetric');
    wave2d_diff_vec(cov_model0,cov_model,m_inds,stpars,0);
    wave2d_diff_vec(cov_imetric0,cov_imetric,m_inds,stpars,0);
    
    % gradient
    %gradient_model(:) = m0(:) / cov_imetric(:)
    gradient_tot           = zeros(nmod,1);
    gradient_model         = (m0-mprior) ./ cov_imetric;
    %gradient_model         = (m0-mprior) ./ cov_model;
    gradient_data          = zeros(nmod,1);
    gradient_data(indB)    = summed_ker .* da_local_vec;
    gradient_data(indTXY)  = source_gradient;
    % zero out parts
    if INV_SOURCE_T==0
        gradient_data(indT) = 0;
        gradient_model(indT) = 0;        
    end
    if INV_SOURCE_X==0
        gradient_data(indX) = 0;
        gradient_model(indX) = 0;  
        gradient_data(indY) = 0;
        gradient_model(indY) = 0; 
    end
    if INV_STRUCT_BETA==0
        gradient_data(indB) = 0;
        gradient_model(indB) = 0;
    end
    gradient_tot = gradient_data + gradient_model;
    
    % check actual norms
    % --> gradient_norm_all.dat, gradient_norm_data_all.dat, gradient_norm_model all.dat
    %wave2d_gnorm_sq(gradient_tot,cov_imetric,m_inds);
    %wave2d_gnorm_sq(gradient_data,cov_imetric,m_inds);
    %wave2d_gnorm_sq(gradient_model,cov_imetric,m_inds);
    wave2d_norm_sq(gradient_tot,mprior,cov_imetric,m_inds,ones(4,1),0);
    wave2d_norm_sq(gradient_data,mprior,cov_imetric,m_inds,ones(4,1),0);
    wave2d_norm_sq(gradient_model,mprior,cov_imetric,m_inds,ones(4,1),0);
    
    % gradient -- TO CHECK
    %write(19,'(3e20.10)') gradient(i), gradient_data(i), gradient_model(i)
    [gradient_tot0,gradient_data0,gradient_model0] = textread([idir1 'gradient_vec.dat'],'%f%f%f');
    % check
    disp('comparing gradient, gradient_model, and gradient_data');
    wave2d_diff_vec(gradient_tot0,gradient_tot,m_inds,stpars,0);
    wave2d_diff_vec(gradient_data0,gradient_data,m_inds,stpars,0);
    wave2d_diff_vec(gradient_model0,gradient_model,m_inds,stpars,0);
     
    % check chi model values
    model_norm_parts = zeros(4,1);
    model_norm_parts(1) = sum( (m0(indB)-mprior(indB)).^2 ./ cov_model(indB) * covm_weight_parts(1) );
    model_norm_parts(2) = sum( (m0(indT)-mprior(indT)).^2 ./ cov_model(indT) * covm_weight_parts(2) );
    model_norm_parts(3) = sum( (m0(indX)-mprior(indX)).^2 ./ cov_model(indX) * covm_weight_parts(3) );
    model_norm_parts(4) = sum( (m0(indY)-mprior(indY)).^2 ./ cov_model(indY) * covm_weight_parts(4) );
    model_norm =  sum(model_norm_parts);
    %model_norm = sum( (m0-mprior).^2 ./ cov_model );
    % check
    model_norm0, model_norm, sum(model_norm_parts)
    
    % check chi model values
    chi_k_val = 0.5*(data_norm + model_norm)
    
    %---------------------------------------------------------
    % emulate CG algorithm in wave2d.f90 -- THIS ASSUMES A DIAGONAL COVARIANCE MATRIX
    
    disp('------CG ALGORITHM---------');
    
    gk = gradient_tot;
    mt = zeros(nmod,1);
    if istep <= 1
        beta_val = 0.0;
        p0 = zeros(nmod,1);
    else
        % step for search direction
        % requires current gradient gk AND previous gradient g0
        
        % load p0 and g0, the PREVIOUS gradient vectors (pk and gk)
        % write(19,'(4e16.6)') g0(i), gk(i), p0(i), pk(i)
        %[~,g0,~,p0] = textread([idir2 'cg_grad_vectors.dat'],'%f%f%f%f');
        [g0,p0] = textread([idir2 gfilemp],'%f%f');
        
        beta_val = sum((gk - g0) .* (cov_imetric .*gk) ) / sum(g0 .* (cov_imetric.*g0) );
        if isinf(beta_val), error('beta_val is infinity'); end
    end
    % search direction vector
    pk = -cov_imetric .* gk + beta_val * p0;
    % step for test model
    mu_val = chi_data_stop;
    lam_t_val_bot = sum( gk .* pk );
    lam_t_val = 2.0*(mu_val - chi_k_val) / lam_t_val_bot;
    
    % check
    if 0==1
        disp('checking values in the CG algorithm');
        vals = [mu_val lam_t_val_bot lam_t_val beta_val chi_k_val];
        %write(19,'(5e16.8)') mu_val, lam_t_val_bot, lam_t_val, beta_val, chi_k_val
        vals0 = load([dirbase 'SEM2D_iterate_OUTPUT/run_9150/cg_test_vals.dat']);
        for ii=1:length(vals), vals0(ii), vals(ii), end
        break
    end
    
    if imaketest==1
        % test model
        mt = m0 + lam_t_val*pk;
        
    else
        % load chi for test model
        chi_t_val = load(chitfile);
        
        % a quadratic fit requires at least 5 values
        xx1 = 0.0
        xx2 = lam_t_val
        yy1 = chi_k_val
        yy2 = chi_t_val
        g1  = sum(gk .* pk);
        %g1  = sum(g0 .* pk);
        
        % coefficients of the quadratic polynomial (ax^2 + bx + c)
        Pa = ((yy2 - yy1) - g1*(xx2 - xx1)) / (xx2^2 - xx1^2)
        Pb = g1
        Pc = yy1 - Pa*xx1^2 - Pb*xx1
        
        % get the analytical minimum (the vertex)
        if (Pa ~= 0.0)
            xmin = -Pb / (2.0*Pa)
        else
            error('check the quadratic input polynomial');
        end
        
        % check
        %vals = load([dirbase 'SEM2D_iterate_OUTPUT/run_1201/quad_poly.dat']);
        %for ii=1:length(vals), vals(ii), end
        
        % compute updated model
        lam_0_val = xmin;
        mk = m0 + lam_0_val * pk;
        
    end

%---------------------------------------------------------

elseif any(icg == [2 3 4])
    
    % load Matlab parameterization of the fields (gaussian_2D.m)
    load([dirrand 'matlab_vars']);
    
    % number of structure parameters
    ncell = nx*ny;
    
    % plotting mesh
    Xp = reshape(x,ny,nx);
    Yp = reshape(y,ny,nx);
    
    % compute the GLL gridpoint (local) closest to each nx x ny cell;
    % this will give a quick mapping from kernels to cells
    %[iGLL, dGLL] = wave2d_gll2cell(xg,yg,x,y);
    
    % indexing for cell-based model
    m_indsi = m_inds - nlocal + ncell; m_indsi(1,1) = 1;    
    indBi = m_indsi(1,1):m_indsi(1,2);
    indTi = m_indsi(2,1):m_indsi(2,2);
    indXi = m_indsi(3,1):m_indsi(3,2);
    indYi = m_indsi(4,1):m_indsi(4,2);
    indTXYi = m_indsi(2,1):m_indsi(4,2);
    
    % vector associated with parameterization
    Atot = nx*ny*dx^2;               % NOTE: does not equal xran*yran
    dA = Atot/n;
    dAvec = ones(n,1)*dA;
    Avec = sqrt( dAvec/Atot );
    n*dA - Atot, sum(Avec.^2)         % check
    
    if 0==1
        % plot full covariance matrix
        figure; imagesc(C);
        title('prior model covariance matrix');
        axis equal, axis tight, colorbar

        % plot ordering of cells for the first tenth of the point
        gcut = round(length(x)/10);
        figure; hold on;
        plot(x(1:gcut),y(1:gcut),'.');
        for ii=1:gcut, text(x(ii),y(ii),num2str(ii),'fontsize',4); end
    end
    
    %--------------------------------------
    % covariance matrix
  
    % modified covariance matrix
    %Cmod  = diag(1./Avec) * C * diag(1./Avec); 
    
    nmodi = ncell + nmod_src;
    Cfull = zeros(nmodi,nmodi);
    Cmod = zeros(nmodi,nmodi);
    cov_weight_vec = zeros(nmodi,1);
    
    if 1==1
        % full, unweighted covariance matrix
        Cfull(indBi, indBi) = C;
        Cfull(indTi, indTi) = diag(sigma_ts^2 * ones(length(indTi),1) );
        Cfull(indXi, indXi) = diag(sigma_xs^2 * ones(length(indXi),1) );
        Cfull(indYi, indYi) = diag(sigma_zs^2 * ones(length(indYi),1) );
        %figure; spy(Cfull);

        cov_weight_vec = zeros(nmodi,1);
        cov_weight_vec(indBi) = coverage_str / covg_weight_parts(1) ./ Avec.^2;
        cov_weight_vec(indTi) = dnevent_run * coverage_src / covg_weight_parts(2);
        cov_weight_vec(indXi) = dnevent_run * coverage_src / covg_weight_parts(3);
        cov_weight_vec(indYi) = dnevent_run * coverage_src / covg_weight_parts(4);
        
        % KEY: weighting for full covariance matrix
        Cmod = diag(sqrt(cov_weight_vec)) * Cfull * diag(sqrt(cov_weight_vec));

    else
        % diagonal, unweighted covariance matrix
        Cfull(indBi, indBi) = diag(diag(C) ./ Avec.^2 * coverage_str );
        Cfull(indTi, indTi) = diag(sigma_ts^2 * dnevent_run * coverage_src * ones(length(indTi),1) );
        Cfull(indXi, indXi) = diag(sigma_xs^2 * dnevent_run * coverage_src * ones(length(indXi),1) );
        Cfull(indYi, indYi) = diag(sigma_zs^2 * dnevent_run * coverage_src * ones(length(indYi),1) );
        
        % test with a diagonal matrix
        Cmod(indBi, indBi) = Cfull(indBi, indBi) / covg_weight_parts(1);
        Cmod(indTi, indTi) = Cfull(indTi, indTi) / covg_weight_parts(2);
        Cmod(indXi, indXi) = Cfull(indXi, indXi) / covg_weight_parts(3);
        Cmod(indYi, indYi) = Cfull(indYi, indYi) / covg_weight_parts(4);
    end
    
    % equivalent representations
    %cdiagmod = diag( diag(sqrt(cov_weight_vec)) * diag(diag(Cfull)) * diag(sqrt(cov_weight_vec)) );
    %cmoddiag = diag(Cmod);
    
    % C-inverse for the model norm computations
    Cmodinv = diag(1./diag(Cmod));
%     if icg==2
%         %cinv = 1./diag(Cfull);
%         %cmodinv = 1./diag(Cmod);
%         %Cinv = diag(1./diag(Cfull));
%         Cmodinv = diag(1./diag(Cmod));
%     else
%         %Cinv = inv(Cfull);
%         Cmodinv = inv(Cmod);
%     end

    if icg==2
        Crun = diag(diag(Cmod));
    else   % icg = 3,4
        Crun = Cmod;
    end
        
    %--------------------------------------
    % models
    
    m0i      = wave2d_m_gll2cell(m0,xg,yg,Xp,Yp,nx,ny,nmod_src,1);
    mpriori  = wave2d_m_gll2cell(mprior,xg,yg,Xp,Yp,nx,ny,nmod_src,1);
    mtargeti = wave2d_m_gll2cell(mtarget,xg,yg,Xp,Yp,nx,ny,nmod_src,1);
    
    %--------------------------------------
    % gradient
    
    % interpolate kernels onto cells
    kernels_alli = zeros(ncell,nevent_run);
    for ii=1:nevent_run 
        % kernels_all is computed above
        Kbeta0 = kernels_all(:,ii);
        %Kbeta = Kbeta0(iGLL);  % crude interpolation
        Kbeta = wave2d_gll2cell(xg,yg,Kbeta0,Xp,Yp);
        kernels_alli(:,ii) = Kbeta(:);

        if iplotker==1
            % plot ghat and g = Cm*ghat
            figure; nr=2; nc=1;
            Zp = Kbeta;
            %Zp = reshape( Kbeta, ny, nx);   % gradient
            subplot(nr,nc,1); imagesc(Zp); set(gca,'ydir','normal'); axis equal, axis tight;
            caxis([-1 1]*max(abs(Kbeta(:)))*0.8); colormap(seis); colorbar;
            title(sprintf('kernel %i out of %i',ii,nevent_run));
            
            Zp = reshape( C * Kbeta(:), ny, nx);  % steepest ascent vector
            subplot(nr,nc,2); imagesc(Zp); set(gca,'ydir','normal'); axis equal, axis tight;
            caxis([-1 1]*max(abs(C * Kbeta(:)))*0.8); colormap(seis); colorbar;
            title(sprintf('kernel %i out of %i',ii,nevent_run));
            orient tall, wysiwyg;
            
            if 0==1
                %print(gcf,'-depsc',sprintf('kernel_%2.2i',ii) );
                %print(gcf,'-dpsc',sprintf('kernel_%2.2i',ii) );
                print(gcf,'-dpdf',sprintf('kernel_%2.2i',ii) );
                pause(1);
            end
        end
    end
    summed_keri = sum(kernels_alli,2);
    
    gradienti_data = zeros(nmodi,1);
    gradienti_model = zeros(nmodi,1);
    gradienti_tot = zeros(nmodi,1);
    gradienti_data(indBi) = summed_keri .* dAvec;
    gradienti_data(indTXYi) = source_gradient;
    gradienti_model = Cmodinv * (m0i - mpriori);
    %gradienti_model = Cinv * (m0i - mpriori);
    % zero out parts
    if INV_SOURCE_T==0
        gradienti_data(indTi) = 0;
        gradienti_model(indTi) = 0;        
    end
    if INV_SOURCE_X==0
        gradienti_data(indXi) = 0;
        gradienti_model(indXi) = 0;  
        gradienti_data(indYi) = 0;
        gradienti_model(indYi) = 0; 
    end
    if INV_STRUCT_BETA==0
        gradienti_data(indBi) = 0;
        gradienti_model(indBi) = 0;
    end
    gradienti_tot = gradienti_data + gradienti_model;
    
    % gradient -- TO CHECK
    % NOTE: due to the new discretization of the kernels, we expect the
    %       norm of the structure gradient to be slightly different from before
    %wave2d_gnorm_sq(gradienti_tot,Crun,m_indsi);
    %wave2d_gnorm_sq(gradienti_data,Crun,m_indsi);
    %wave2d_gnorm_sq(gradienti_model,Crun,m_indsi);
    wave2d_norm_sq(gradienti_tot,mpriori,Crun,m_indsi,ones(4,1),0);
    wave2d_norm_sq(gradienti_data,mpriori,Crun,m_indsi,ones(4,1),0);
    wave2d_norm_sq(gradienti_model,mpriori,Crun,m_indsi,ones(4,1),0);
    
    % check parts of model norm
    model_norm_partsi = wave2d_norm_sq(m0i,mpriori,Crun,m_indsi,covm_weight_parts,1);
    model_normi = sum(model_norm_partsi);
    % check
    model_norm0, model_normi % sum( m0i(1:ncell).^2 .* cinv(1:ncell) )
    
    % check chi model values
    chi_k_val = 0.5*(data_norm + model_normi)
    chi_k_val = chi0        % just use the real chi value
    
    %---------------------------------------------------------
    
    if icg ~=4 
        % emulate CG algorithm in wave2d.f90 -- THIS ASSUMES A FULL COVARIANCE MATRIX

        % compute new model vector (mt or mk)
        gfile = [idir2 gfilemp];
        mu_val = chi_data_stop;
        gk = gradienti_tot;
        [M,pk] = wave2d_cg(m0i,gk,Crun,chi_k_val,mu_val,istep,imaketest,gfile,chitfile);

        if imaketest==1
            mti = M;
        else
            mki = M;
        end
    
    else
        % blocks from Cm
        Cstr = Cmod(indBi,indBi);
        csrc = diag(Cmod(indTXYi,indTXYi));   % diagonal
        
        % gradient for structure
        Gstr = zeros(nevent_run,ncell);
        wvec = sqrt( dnorm2 );
        Gstr = -kernels_alli' .* repmat(dAvec',nevent_run,1) .* repmat(1./wvec, 1, ncell);
        
        % gradient for source
        Gsrc = zeros(nevent_run,nmod_src);
        % QUESTION: SHOULD THIS BE THE FULL GRADIENT OR JUST THE DATA TERM?
        %Gsrc = [diag(gradienti_tot(indTi)) diag(gradienti_tot(indXi)) diag(gradienti_tot(indYi))];
        
        % partial derivatives matrix (nsrc x nmodi)
        G = [Gstr Gsrc];
        
        % modified data vector
        dmod = sqrt(dnorm2) + G*(m0i-mpriori);
        
        % compute subspace Hessian
        [Hstr,Hsrc] = wave2d_compute_hessian(Gstr,Gsrc,Cstr,csrc);
        H = Hstr + Hsrc + eye(nsrc,nsrc);
        wave2d_plot_hessian(Hstr,Hsrc,H);
        
        % write out files for Hessian
        %if iwrite==1
        %    nprint = 20;
        %    wave2d_write_hessian(H,eids,dir2,nprint);
        %end
        
        % examining H\dnorm2 operation with TSVD
        %mu_all = wave2d_tsvd_hessian(H,dnorm2);
        %mu_all = wave2d_tsvd_hessian(H,sqrt(dnorm2));
        mu_all = wave2d_tsvd_hessian(H,dmod);
        if 0==1
            lambda = 1;
            ndiff = zeros(nsrc,1);
            for ii=1:nsrc
            %for ii=1:5    
                mu = mu_all(:,ii);
                mplot = lambda * Cstr * Gstr' * mu;
                pmax = max(abs(mplot));
                pmax = 0.1;
                
                ndiff(ii) = norm( mplot - (mtargeti(indBi) - m0i(indBi)) );
                
                if 1==1
                    Zp = griddata(x,y,mplot,Xp,Yp,'nearest');
                    figure; nr=2; nc=1;
                    subplot(nr,nc,1); imagesc(Zp); set(gca,'ydir','normal');
                    title(sprintf('dm -- TSVD truncated at index %i / %i, norm = %.2e',ii,nsrc,ndiff(ii)));
                    caxis([-1 1]*pmax); colormap(seis); colorbar;
                    axis equal, axis tight

                    subplot(nr,nc,2); hold on;
                    if ii > 1, plot(mu_all(:,ii-1),'ro--','markersize',8,'linewidth',1); end
                    if ii < nsrc, plot(mu_all(:,ii+1),'ko--','markersize',8,'linewidth',1); end
                    plot(mu,'b.-','markersize',20,'linewidth',2);
                    axis([0 nsrc+1 max(abs(mu))*[-1 1]]);
                    grid on; xlabel('event index'); ylabel('elements of mu vector');
                    title('\mu_{-1} (RED)  -->  \mu_0 (BLUE)  -->  \mu_{+1} (BLACK)');
                    orient tall, wysiwyg

                    print(gcf,'-dpdf',sprintf('%s_m%s_dm_tsvd_%2.2i',stirun0,stim0,ii) );
                end
            end
            
            % difference between target model and new model
            figure; plot(ndiff,'.-','markersize',20,'linewidth',2); grid on;
            xlabel('TSVD truncation index'); ylabel('norm( dm0 - (mtar-m0) )');
            print(gcf,'-dpdf',sprintf('%s_m%s_dm_ndiff',stirun0,stim0) );
        end
        
        % initialize to no update
        dm_str = zeros(nmod_str,1);
        dm_src = zeros(nmod_src,1);
        dm_str1 = zeros(nmod_str,1); dm_str2 = zeros(nmod_str,1);
        dm_src1 = zeros(nmod_src,1); dm_src2 = zeros(nmod_src,1);
        
        lambda = 1;
        itest = 2;
        switch itest
            case 1
                % KEY: compute model update
                mu = H \ sqrt(dnorm2);    % why not dmod?
                %mu = mu_all(:,15);
                
                % compute the STRUCTURE model update
                if INV_STRUCT_BETA == 1
                    dm_str = Cstr * Gstr' * mu;
                end
                % compute the SOURCE model update
                if and(INV_SOURCE_T==1, INV_SOURCE_X==1)
                    dm_src = csrc .* (Gsrc' * mu);
                end
                dm = [dm_str ; dm_src];
                mki = m0i + dm;
                
            case 2  % actual iterative algorithm
                mu = H \ dmod;
                if INV_STRUCT_BETA == 1
                    dm_str = Cstr * Gstr' * mu;
                end
                if and(INV_SOURCE_T==1, INV_SOURCE_X==1)
                    dm_src = csrc .* (Gsrc' * mu);
                end
                dm = [dm_str ; dm_src];
                mki = mpriori + dm;
         
            case 3  % actual iterative algorithm
                if INV_STRUCT_BETA == 1
                    mu1 = H \ sqrt(dnorm2);
                    mu2 = H \ (Gstr*(m0i(indBi)-mpriori(indBi)));
                    dm_str1 = Cstr * Gstr' * mu1;
                    dm_str2 = Cstr * Gstr' * mu2;
                end
                dm1 = [dm_str1 ; dm_src1];
                dm2 = [dm_str2 ; dm_src2];
                dm = dm1 + dm2;
                mki = mpriori + dm;         % mpriori or m0i
                figure; plot( dm1, dm2 ,'.')
                
            case 4
                if INV_STRUCT_BETA == 1
                    dm_str1 = (Gstr'*Gstr + Cmodinv(indBi,indBi) ) \ (Gstr'*sqrt(dnorm2));
                    dm_str2 = -(diag(diag(Cstr))*Gstr'*Gstr + eye(ncell,ncell) ) \ (m0i(indBi)-mpriori(indBi));
                end
                dm1 = [dm_str1 ; dm_src1];
                dm2 = [dm_str2 ; dm_src2];
                dm = dm1 + dm2;
                mki = m0i + dm;
                figure; plot( dm1, dm2 ,'.')
        end
        
        % linear combination of kernels
        figure; plot(mu,'.-'); xlabel('event index');
        title('mu vector for source subspace approach');
        
        figure; nr=3; nc=3; cmax = 0.1;
        for ii=1:6
            switch ii
                case 1, mplot = mpriori(1:ncell); smt = 'mprior';
                case 2, mplot = mtargeti(1:ncell); smt = 'mtarget'; 
                case 3, mplot = m0i(1:ncell); smt = 'm0';
                case 4, mplot = mki(1:ncell); smt = 'mnew';
                    if itest==3, smt = 'mnew = mpriori + dm'; else smt = 'mnew = m0 + dm'; end
                case 5, mplot = mki(1:ncell) - m0i(1:ncell); smt = 'mnew - m0'; 
                case 6, mplot = mtargeti(1:ncell) - m0i(1:ncell); smt = 'mtarget - m0'; 
            end
            Zp = griddata(x,y,mplot,Xp,Yp,'nearest');
            subplot(nr,nc,ii);
            imagesc(Zp); set(gca,'ydir','normal'); title(smt)
            caxis([-1 1]*cmax); colormap(seis); axis equal, axis tight
        end
        
        if or(itest==3, itest==4)
            for ii=7:9
                switch ii
                    case 7, mplot = dm(1:ncell); smt = 'dm';
                    case 8, mplot = dm1(1:ncell); smt = 'dm(term 1)';
                    case 9, mplot = dm2(1:ncell); smt = 'dm(term 2)';
                end
                Zp = griddata(x,y,mplot,Xp,Yp,'nearest');
                subplot(nr,nc,ii);
                imagesc(Zp); set(gca,'ydir','normal'); title(smt)
                caxis([-1 1]*cmax); colormap(seis); axis equal, axis tight
            end
            orient tall, wysiwyg
        else
            % uniform sum to make 'misfit kernel' -- this is the pattern for the CG algorithm
            mplot = Cstr * Gstr' * ones(nsrc,1);
            pmax = max(abs(mplot));
            Zp = griddata(x,y,mplot,Xp,Yp,'nearest');
            subplot(nr,nc,7); imagesc(Zp); set(gca,'ydir','normal'); title('dm (CG)');
            caxis([-1 1]*pmax); colormap(seis); axis equal, axis tight
            orient tall, wysiwyg            
        end
        
%         disp(' norm, min, max of structure update:');
%         norm_dm_str = dm_str' * diag(1./diag(Cstr)) * dm_str
%         min_dm_str = min( dm_str )
%         max_dm_str = max( dm_str )
%         disp(' norm, min, max of source update:');
%         norm_dm_src = sum( dm_src .* dm_src ./ csrc )
%         min_dm_src = min( dm_src )
%         max_dm_src = max( dm_src )
    end
    
 
else
    error('invalid icg value');
    
end

%========================================================================
% WRITE NEW MODEL TO FILE

if icg ~= 4     % CG method
    % write current gradient and search direction to file
    if iwrite==1
        wave2d_write_grad([idir1 gfilem0],gk,pk);  % note idir1
    end

    if icg==1
        if imaketest==1, M = mt; else M = mk; end
        M0 = m0; Mtar = mtarget;
        x0 = xg; y0 = yg;
    else
        if imaketest==1, M = mti; else M = mki; end
        M0 = m0i; Mtar = mtargeti;
        x0 = x; y0 = y;
        m_inds = m_indsi;
    end

else            % source subspace method
    M = mki;
    M0 = m0i; Mtar = mtargeti;
    x0 = x; y0 = y; 
    m_inds = m_indsi;
end

% convert to physical parameters (only differes for structure parameters)
% --> current, new, and target models
M0_vec = wave2d_m2mvec(M0,m_inds,beta0);
M_vec = wave2d_m2mvec(M,m_inds,beta0);
Mtar_vec = wave2d_m2mvec(Mtar,m_inds,beta0);

% split into components
[m0_B, m0_ts, m0_xs, m0_ys] = wave2d_splitm(M0,m_inds);
[m_B, m_ts, m_xs, m_ys] = wave2d_splitm(M,m_inds);
[mtar_B, mtar_ts, mtar_xs, mtar_ys] = wave2d_splitm(Mtar,m_inds);

[m0_vec_B, m0_vec_ts, m0_vec_xs, m0_vec_ys] = wave2d_splitm(M0_vec,m_inds);
[m_vec_B, m_vec_ts, m_vec_xs, m_vec_ys] = wave2d_splitm(M_vec,m_inds);
[mtar_vec_B, mtar_vec_ts, mtar_vec_xs, mtar_vec_ys] = wave2d_splitm(Mtar_vec,m_inds);

if iplotmod==1
    % change w.r.t. previous model
    %Zp = griddata(x0,y0,log(m_vec_B ./ m0_vec_B),Xp,Yp,'nearest');
    %figure; imagesc(Zp); set(gca,'ydir','normal');
    %colormap(seis); caxis([-1 1]*0.1); colorbar;
    
    figure; nr=3; nc=2; cmax = 0.1;
    if imaketest==1
        tlabs = {'ln(m0/beta0)','mtar - m0','ln(mt/beta0)','mtar - mt','ln(mtar/beta0)','mt - m0'};
    else
        tlabs = {'ln(m0/beta0)','mtar - m0','ln(mk/beta0)','mtar - mk','ln(mtar/beta0)','mk - m0'};
    end
    for kk = 1:6
        switch kk
            case 1, mplot = log(m0_vec_B/beta0);
            case 2, mplot = log(mtar_vec_B./m0_vec_B);
            case 3, mplot = log(m_vec_B/beta0);
            case 4, mplot = log(mtar_vec_B./m_vec_B);
            case 5, mplot = log(mtar_vec_B/beta0);   
            case 6, mplot = log(m_vec_B./m0_vec_B);   
        end
        
        Zp = griddata(x0,y0,mplot,Xp,Yp,'nearest');
        %mnorm = sum( mplot.^2 ./ cov_beta' );

        subplot(nr,nc,kk); hold on;
        imagesc(Zp); set(gca,'ydir','normal');
        caxis([-1 1]*cmax); colormap(seis); axis equal, axis tight
        title(tlabs{kk});
        %plot(rlon,rlat,'k.','markersize',16);
        %plot(slon,slat,'p','markersize',18,'markeredgecolor','k','linewidth',2);
        %title([tlabs{kk} sprintf(', norm = %.2f',mnorm)]);
    end
    orient tall, wysiwyg
end

if iwrite==1
    % interpolate onto SEM mesh
    if icg ~= 1
        m_B = wave2d_cell2gll(xg,yg,Xp,Yp,reshape(m_B,ny,nx));
    end
    
    % write structure model -- ONLY MU CHANGES (only m_B appears)
    m_str_beta_new  = beta0 * exp( m_B );
    m_str_rho_new   = rho0 * ones(nlocal,1);
    m_str_kappa_new = kappa0 * ones(nlocal,1);
    m_str_mu_new    = m_str_rho_new .* m_str_beta_new.^2;
    wave2d_write_str([odir 'structure_syn_m' stimo '.dat'],xg,yg,...
        m_str_kappa_new,m_str_mu_new,m_str_rho_new,m_B);
    
    % check
    %mcheck = mt_vec0(indB);
    %norm(mcheck), norm(m_str_beta_new), norm(mcheck - m_str_beta_new)/norm(m_str_beta_new)
    
    % read sources for data
    [m_src_lon,m_src_lat,m_src_ts_dat,m_src_xs_dat,m_src_ys_dat,m_src_ts_d,m_src_xs_d,m_src_ys_d] ...
        = textread([idir1 'src_dat.dat'],'%f%f%f%f%f%f%f%f');
    % mdat - msrc : how far to still go
    m_src_ts_d_new = m_src_ts_dat - m_vec_ts;
    m_src_xs_d_new = m_src_xs_dat - m_vec_xs;
    m_src_ys_d_new = m_src_ys_dat - m_vec_ys;    
    wave2d_write_src([odir 'src_syn_m' stimo '.dat'], m_src_lon, m_src_lat,...
        m_vec_ts, m_vec_xs, m_vec_ys, m_src_ts_d_new, m_src_xs_d_new, m_src_ys_d_new);
    
%     % check all
%     if icg==1
%         mcheck = mt_vec0;
%         for ii=1:4
%             inds = m_inds(ii,1):m_inds(ii,2);
%             disp('------');
%             disp(sprintf('%.3e %.3e %.3e %.3e',...
%             norm(mcheck(inds)), norm(M_vec(inds)),...
%             norm(mcheck(inds) - M_vec(inds))/norm(mcheck(inds)) ));
%         end
%     end
    
end

%=============================================================
