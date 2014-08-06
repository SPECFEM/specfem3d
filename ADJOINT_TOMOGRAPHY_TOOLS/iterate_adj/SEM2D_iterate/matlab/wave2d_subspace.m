%
% wave2d_subspace.m
% Carl Tape, 26-Jan-2010
%
% This program implements the subspace method inversion notes developed by
% Malcolm Sambridge, Jeroen Tromp, and Carl Tape.
%
% This assumes a diagonal model covariance matrix and requires
% regularization -- several choices are explored here.
% 
% calls xxx
% called by xxx
%

clc, clear
%close all
format short, format compact
%warning off

% add path to additional matlab scripts
path([pwd '/matlab_scripts'],path);

colors;

ax1 = [-121 -114 31 37];
stfm = '%4.4i';

%----------------------------------------------
% USER INPUT

dir_run = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_OUTPUT/';

TSVD = 1;       % truncated singular value decomposition
inew = 0;       % subspace AND parameter class inversions

%----------------------------------------------

irun0 = input(' Enter irun0 : ');
iread = input(' Enter 1 to read in new event kernels (0 otherwise) : ');
hmod = input(' Enter next model number (hmod = 2,4,6,8,...) : ');
INV_STRUCT = input(' Enter 1 to invert for STRUCTURE : ');
INV_SOURCE = input(' Enter 1 to invert for SOURCE : ');
iwrite = input(' Enter 1 to write out files : ');
if iwrite==1
    COMPUTE_KERNELS = input(' Enter 1 to compute the event kernels for this model : ');
end
if INV_STRUCT==1
    ifigure = input(' Enter 1 to plot figures : ');
else
    ifigure = 0;
end

%iwrite = 0;
%ifigure = 1;
%hmod = 1;   % 1, 2, 3, ...
%irun0 = 7150;     % 6080, 6180, 7000, 7050

NLOCAL = 40000;     % assume 1 source parameter only for now (beta)
NPARM_SOURCE = 3;
NPARM_STRUCT = 1;

nsrc = 25;
nrec = 132;
nmeas = nsrc * nrec;

nmod_str = NLOCAL;
nmod_src = nsrc * NPARM_SOURCE;
nmod = nmod_str + nmod_src;

%SIGMA_DT = 0.20;
%sval_cut = SIGMA_DT * sqrt(nmeas);
sval_cut = 5;

% event IDs
for ii=1:nsrc, eids(ii) = cellstr(sprintf(stfm,ii)); end

%----------------

stm = sprintf(stfm,hmod);

% directories
dir0 = [dir_run 'run_' sprintf(stfm,irun0) '/'];
if ~exist(dir0), error([dir0 ' does not exist']); end
dir0r = [dir0 'READ_IN/'];
if ~exist(dir0r), mkdir(dir0r); end
dir2 = [dir0r 'model_m' stm '/'];
if ~exist(dir2), mkdir(dir2); end
dir2lab = ['run-' sprintf(stfm,irun0) ', READ-IN -- toward model m' stm ];

if hmod == 2
    dir1 = dir0;
else
    dir1 = [dir0 'READ_IN/model_m' sprintf(stfm,hmod-2) '/'];
end

% indexing for source parameters
itemp = NPARM_SOURCE*[1:25]';
indmat2 = [ [1 ; 1+itemp(1:end-1)] itemp ];

% indexing for model vector
m_inds = load([dir1 'm_inds.dat']);
inds_B  = [m_inds(1,1) : m_inds(1,2)];
inds_ts = [m_inds(2,1) : m_inds(2,2)];
inds_xs = [m_inds(3,1) : m_inds(3,2)];
inds_ys = [m_inds(4,1) : m_inds(4,2)];
inds_src = [m_inds(2,1) : m_inds(4,2)];

% load the source parameters
src_syn = load([dir1 'src_syn.dat']);
slon = src_syn(:,1);
slat = src_syn(:,2);
src_ts = src_syn(:,3);
src_xs = src_syn(:,4);
src_ys = src_syn(:,5);

% load the source parameters
src_dat = load([dir1 'src_dat.dat']);
src_ts_dat = src_dat(:,3);
src_xs_dat = src_dat(:,4);
src_ys_dat = src_dat(:,5);

% load the receivers
rec_lonlat = load([dir1 'recs_lonlat.dat']);
rlon = rec_lonlat(:,1);
rlat = rec_lonlat(:,2);

% load the covariance matrix
% --> cov_imetric(NLOCAL+1 : nmod_str) = ( sigma_beta  )**2 / da_local_vec(:) * AREA
%cov_model = load('/net/denali/scratch1/carltape/OUTPUT_2/run_9100/cov_imetric_diagonal.dat');
[cov_model, cov_imetric] = textread([dir1 'cov_model_diagonal.dat'],'%f%f');
cov_beta  = cov_model(inds_B)';
cov_src   = cov_model(inds_src)';
clear cov_model_all

% load the reference values
% alpha0, beta0, rho0, bulk0, kappa0, mu0
vall = load([dir1 'reference_values.dat']);
alpha0 = vall(1);
beta0  = vall(2);
rho0   = vall(3);
bulk0  = vall(4);
kappa0 = vall(5);
mu0    = vall(6);

% load the model vector
mtemp = load([dir1 'cg_model_vectors.dat']);
m_all = mtemp(:,1);
m_src   = mtemp(inds_src)';

% prior, initial, and target models
[mprior, m0_initial, mtarget] = textread([dir0 'prior_initial_target_models.dat'],'%f%f%f');

% load the structure files
mtemp = load([dir1 'structure_syn.dat']);
m_str_lon   = mtemp(:,1);
m_str_lat   = mtemp(:,2);
m_str_kappa = mtemp(:,3);
m_str_mu    = mtemp(:,4);
m_str_rho   = mtemp(:,5);
m_str_B     = mtemp(:,6);

% load the structure files
mtemp = load([dir1 'structure_dat.dat']);
m_dat_str_lon   = mtemp(:,1);
m_dat_str_lat   = mtemp(:,2);
m_dat_str_kappa = mtemp(:,3);
m_dat_str_mu    = mtemp(:,4);
m_dat_str_rho   = mtemp(:,5);
m_dat_str_B     = mtemp(:,6);

% load source files (synthetics and data)
[m_src_lon,m_src_lat,m_src_ts,m_src_xs,m_src_ys,m_src_ts_d,m_src_xs_d,m_src_ys_d] ...
    = textread([dir1 'src_syn.dat'],'%f%f%f%f%f%f%f%f');
[junk1,junk2,m_src_ts_dat,m_src_xs_dat,m_src_ys_dat,junk3,junk4,junk5] ...
    = textread([dir1 'src_dat.dat'],'%f%f%f%f%f%f%f%f');

% load the gradient
gtemp = load([dir1 'gradient_vec.dat']);
gradient = gtemp(:,1);
grad_beta = gradient(inds_B)';
grad_src  = gradient(inds_src)';
disp(' check norms of structure and source GRADIENTS :');
norm_grad_str = dot( grad_beta, cov_beta.*grad_beta )
norm_grad_src = dot(  grad_src,  cov_src.*grad_src )
disp(' gradient balance for the CG inversion :');
norm_grad_src / norm_grad_str

% load the source gradient and partition into matrix
%grad_src = load([dir1 'source_gradient.dat']);
% grad_src_mat = zeros(nsrc,NPARM_SOURCE);
% for isrc = 1:nsrc
%     grad_src_mat(isrc,:) = grad_src(indmat2(isrc,1) : indmat2(isrc,2));
% end

% check the misfit function value
chi = load([dir1 'chi.dat']);
data_norm = load([dir1 'data_norm.dat']);
model_norm = load([dir1 'model_norm.dat']);

%====================================================

% load the data covariance matrix
cov_data = load([dir1 'cov_data_diagonal.dat']);

% load all the measurements and partition into matrix
meas_all = load([dir1 'measure_vec.dat']);
dT_all = meas_all(:,1);     % WITH ERRORS ADDED

% vector of number of measurement windows per event
Ns_vec = nrec * ones(nsrc,1);   % KEY: same measurements per event

% KEY: compute projected data vector
[dnorm, indmat] = wave2d_dsub(dT_all,cov_data,Ns_vec);

% compute the weights (SIGN OR NOT?)
%ws = zeros(nsrc,1);
%ws = 1 ./ sqrt( sum( dT_norm_mat, 2) );
%if sum(isinf(ws)) > 0, error(' For at least one source, there is perfect fit.'); end
% compute the new data vector
%dnorm = 1 ./ ws;

disp('  '); disp(' CHECKING VARIOUS NORMS:');
disp('   model norm:');
model_norm, sum( (m_all-mprior).^2 ./ cov_model )
disp('   data norm:');
data_norm, sum(dnorm)
disp('   misfit function value:');
chi, 0.5*( sum( dnorm ) + sum( (m_all-mprior).^2 ./ cov_model ))

%--------------------------------------------------------------------------
% ASSEMBLING THE GRADIENTS

% initialize
Gstr = zeros(nsrc,NLOCAL);
Gsrc = zeros(nsrc,nmod_src);

if INV_STRUCT == 1

    % load the jacobian for constructing the "event gradient" from the event kernel
    lmesh_all = load([dir1 'local_mesh.dat']);
    Ai = lmesh_all(:,9)';

    % load the event kernels
    %iread = 1;
    efile = 'wave2d_kernel';
    if iread==1
        ismooth = input(' Enter 1 to read smoothed event kernels (0 otherwise) : ');
        disp('reading in the event kernels...');
        Kall = zeros(nsrc,NLOCAL);
        for isrc = 1:nsrc
            isrc
            dirK = [dir1 'event_' sprintf('%3.3i',isrc) '/'];
            kernel = load([dirK 'kernel_basis_smooth']); Kbeta = kernel(:,3)';
            if ismooth==0, Kbeta = kernel(:,4); else Kbeta = kernel(:,3); end
            %if ismooth == 1
            %    kernel = load([dirK 'kernel_basis_smooth']); Kbeta = kernel(:,3)';
            %else
            %    kernel = load([dirK 'kernel_basis']); Kbeta = kernel(:,7)';
            %end
            lon = kernel(:,1);
            lat = kernel(:,2);
            Kall(isrc,:) = Kbeta;
        end
        save(efile,'Kall','lon','lat');
        %break
    else
        load(efile);
    end

    % construct G (plot event kernels too)
    %nsrc = 5;
    for isrc = 1:nsrc
        Kbeta = Kall(isrc,:);
        %Gstr(isrc,:) = -ws(isrc) * Kbeta .* Ai;   % SIGN OR NOT?
        Gstr(isrc,:) = -Kbeta .* Ai;

        if 0==1
            [X,Y,Z] = griddataXB(lon,lat,Kbeta,100,'nearest');
            figure; cmax = 1e-7; hold on;
            pcolor(X,Y,Z); shading interp;
            caxis([-1 1]*cmax); colormap(seis);
            plot(rlon,rlat,'k.','markersize',16)
            for irec = 1:nrec
               text(rlon(irec),rlat(irec),sprintf('%6.1f',dT_mat(isrc,irec)),'fontsize',12);
            end
            plot(slon(isrc),slat(isrc),'p','markersize',24,...
                'markerfacecolor','w','markeredgecolor','k','linewidth',2);
        end
    end

%     % Hessian for structure parameters (from event kernels)
%     disp('constructing the Hessian...');
%     Hstr = zeros(nsrc,nsrc);
%     for ii = 1:nsrc
%         for jj = 1:nsrc
%             Hstr(ii,jj) = dot(Gstr(ii,:), cov_beta.*Gstr(jj,:));
%         end
%     end

end

% projected gradient for source parameters 
%Hsrc = zeros(nsrc,nsrc);
%Hsrc_vec = zeros(nsrc,1);
if INV_SOURCE == 1
    
    % QUESTION: SHOULD THIS BE THE FULL GRADIENT OR JUST THE DATA TERM?
    Gsrc = [diag(gradient(inds_ts)) diag(gradient(inds_xs)) diag(gradient(inds_ys))];
    
%     %for isrc = 1:nsrc
%     %    inds = indmat2(isrc,1) : indmat2(isrc,2);
%     %    %Gsrc(isrc,inds) = -ws(isrc) * grad_src(inds);
%     %    Gsrc(isrc,inds) = grad_src(inds);
%     %end
%     Hsrc = Gsrc * diag(cov_src) * transpose(Gsrc);
%     
%     %for isrc = 1:nsrc
%     %    inds = indmat2(isrc,1) : indmat2(isrc,2);
%     %    Hsrc_vec(isrc) = ws(isrc)^2 * dot(grad_src(inds), cov_src(inds).*grad_src(inds) );
%     %end
%     %Hsrc = diag(Hsrc_vec);

end

%--------------------------------------------------------------------------
% COMPUTE THE HESSIAN IN THE SUBSPACE OF SOURCES

[Hstr,Hsrc] = wave2d_compute_hessian(Gstr,Gsrc,cov_beta,cov_src);

% overall Hessian
% NOTE: identity matrix term
H = Hstr + Hsrc + eye(nsrc,nsrc);

% write out files for Hessian
if iwrite==1
    nprint = 20;
    wave2d_write_hessian(H,eids,dir2,nprint);
end

% construct projection matrix (and assign gradient)
if INV_STRUCT == 1
    G = Gstr;
    Pstr = Gstr';
    P = Pstr;
end
if INV_SOURCE == 1
    G = Gsrc;
    Psrc = zeros(nmod_src,nsrc*2);      % 2: origin time and location
    Psrc(1:nsrc,1:nsrc) = Gsrc(1:nsrc,1:nsrc);
    Psrc(nsrc+1:2*nsrc,nsrc+1:2*nsrc) = Gsrc(1:nsrc,nsrc+1:2*nsrc);
    Psrc(2*nsrc+1:3*nsrc,nsrc+1:2*nsrc) = Gsrc(1:nsrc,2*nsrc+1:3*nsrc);
    P = Psrc;
end
if and(INV_STRUCT == 1, INV_SOURCE == 1)
    G = [Gstr Gsrc];
    P = zeros(nmod,3*nsrc);
    P(1:nmod_str,1:nsrc) = Pstr;
    P(nmod_str+1:nmod,nsrc+1:3*nsrc) = Psrc;
end

% construct matrices for joint inversions
if inew == 1
    GG = zeros(nsrc,NCLASS*nsrc);
    CmP = P; for ii=1:nmod, CmP(ii,:) = P(ii,:) * cov_model(nmod_str+ii); end;
    GG = G * CmP;
    Cmi = P' * CmP;

    figure; nr=3; nc=1;
    subplot(nr,nc,1); spy(P(nmod_str-100:nmod_str,:),3); title('bottom 100 rows of P');
    subplot(nr,nc,2); spy(GG,3); title('GG = G Cm P');
    subplot(nr,nc,3); spy(Cmi,3); title('Cmi = P^T Cm P');
    orient tall, wysiwyg
end

disp(' Hessian diagonal contributions from structure and source:');
disp('  structure   source    total    source/structure');
pmat = [diag(Hstr) diag(Hsrc) diag(H) diag(Hsrc)./diag(Hstr) ];
for ii=1:nsrc
    disp(sprintf('%10.2e %10.2e %10.2e %8.2f',pmat(ii,:)))
end
%disp([diag(Hstr) diag(Hsrc) diag(H) diag(Hsrc)./diag(Hstr) ]);
disp(' gradient balance for the Hessian (subspace) inversion (mean of last column) :');
disp(mean( diag(Hsrc)./diag(Hstr) ))

% plot Hessian matrices
wave2d_plot_hessian(Hstr,Hsrc,H);

if INV_SOURCE == 1
   figure; nr=2; nc=1;
   subplot(nr,nc,1); spy(Gsrc); title('Gsrc');
   subplot(nr,nc,2); spy(Gsrc*Gsrc'); title('Gsrc * Gsrc^T');
end

% set of pmax TSVD values to try
iyes = 0;
while iyes == 0
    iyes = input(' Enter 1 if this balance looks okay : ');
end

% check the balance of the gradients -- SAME AS CHECKING THE HESSIAN DIAGONAL
if 0 == 1
    cov_beta0 = cov_beta;
    cov_src0 = cov_src;
    
    kvec = linspace(0.2,6,100);
    for k = 1:length(kvec)
        fac = kvec(k);
        cov_beta = cov_beta0 * fac;
        
        for ii = 1:nsrc
           norm_grad_str(ii) = sum( Gstr(ii,:).^2 .* cov_beta );
           norm_grad_src(ii) = sum( Gsrc(ii,:).^2 .* cov_src );
           norm_grad_tot(ii) = norm_grad_str(ii) + norm_grad_src(ii);
        end
        %disp(' Norms of the gradients and constituent parts:');
        %disp('  structure   source    total    source/structure');
        %disp([norm_grad_str' norm_grad_src' norm_grad_tot' norm_grad_src'./norm_grad_str']);
        %disp(' mean of the last column :');
        disp([fac mean( norm_grad_src'./norm_grad_str' )]);
    end
        
end

%------------------------------------------------------------------------

% truncated singular value decomposition
if TSVD == 1
    
    % analyze the singular values of H
    [U,S,V] = svd(H);
    s = diag(S);
    
    if 0==1
        % analyze the singular values of H
        % (See also tsvd.m)
        [U,S,V] = svd(H);
        s = diag(S);
        p = sum( s > sval_cut );     % KEY: singular value truncation index
        sp = s(1:p);
        Sp = diag(sp);
        Up = U(:,1:p);
        Vp = V(:,1:p);
        whos U S V
        whos Up Sp Vp
        Hp = Up*Sp*Vp';
        Hinv = Vp*inv(Sp)*Up';
        mu = Hinv * dnorm;
        
        if 0==1
            Ncheck = zeros(nsrc,1);
            for p = 1:nsrc
                sp = s(1:p);
                Sp = diag(sp);
                Up = U(:,1:p);
                Vp = V(:,1:p);
                Hp = Up*Sp*Vp';

                Hinv = Vp*inv(Sp)*Up';
                mu = Hinv * dnorm;
                dm_str = transpose(Gstr) * mu .* cov_beta';  % order?
                m_str_B_new = m_str_B + dm_str;
                Ncheck(p) = sum( (dnorm - Gstr*dm_str).^2 );  % matches rss from tsvd.m
            end
        end
        
        % check the norms and the norms of the inverse
        norm(H), norm( U*S*V' ), norm(Hp)
        norm(inv(H)), norm( Vp*inv(Sp)*Up' ), norm(inv(Hp))
        
        figure; hold on; plot(H(:),'b.'); plot(Hp(:),'ro');
    end
    
    %H = H / 0.5;   % 9550
    
    error('replace the following code with wave2d_tsvd_hessian.m');
    
    pinds = [1:nsrc]';
    [mu_all,rss,f_r_ss] = tsvd(dnorm,H,pinds);   % KEY: TSVD
    
    % norms of mu vectors
    mu_norm = zeros(nsrc,1);
    for ip = 1:nsrc
        mu_norm(ip) = norm(mu_all(:,ip));
    end
    
    figure; nr=2; nc=2;
    xlab1 = 'p, singular value index';
    xlab2 = 'p, singular value truncation index';
    ylab1 = 'singular value';
    ylab2 = 'misfit : dot[ d - H*mu(p), d - H*mu(p) ]';
   
    subplot(nr,nc,1); plot(pinds,s,'.','markersize',20);
    grid on; xlabel(xlab1); ylabel(ylab1); title(dir2lab);
    subplot(nr,nc,2); semilogy(pinds,s,'.','markersize',20);
    grid on; xlabel(xlab1); ylabel(ylab1);
    subplot(nr,nc,3); plot(pinds,rss,'.','markersize',20);
    grid on; xlabel(xlab2); ylabel(ylab2);
    subplot(nr,nc,4); semilogy(pinds,rss,'.','markersize',20);
    grid on; xlabel(xlab2); ylabel(ylab2);
    orient tall, wysiwyg
    
    figure; nr=2; nc=2;
    ylab3 = 'norm of mu vector';
   
    subplot(nr,nc,1); semilogy(pinds,s,'.','markersize',20);
    grid on; xlabel(xlab1); ylabel(ylab1);
    subplot(nr,nc,2); semilogy(pinds,rss,'.','markersize',20);
    grid on; xlabel(xlab2); ylabel(ylab2);
    subplot(nr,nc,3); plot(pinds,mu_all,'.');
    grid on; xlabel('source index'); ylabel('elements of mu vectors');
    subplot(nr,nc,4); semilogy(pinds,mu_norm,'.','markersize',20);
    grid on; xlabel(xlab2); ylabel(ylab3);
    orient tall, wysiwyg
else
    Hinv = inv(H);
    mu = Hinv * dnorm;
end

%------------------------------------------------------------------------

% set of pmax TSVD values to try (if TSVD = 1)
if TSVD == 1
    iyes = 0;
    while iyes == 0
        disp(' Enter truncation values for TSVD:');
        pmin = input(' Enter pmin (1): ');
        pmax = input([' Enter pmax (' num2str(nsrc) '): ']);
        pinc = input(' Enter pinc (1): ');
        pmax_vec = [pmin : pinc : pmax]
        iyes = input(' Enter 1 if this is what you want (0 otherwise) : ');
    end
else
    pmax_vec = 1;
end
nump = length(pmax_vec);

source_error_norm_mat = zeros(nump,5);
source_error_norm_mat(:,1) = pmax_vec';

if inew == 1
    Hnew = GG'*GG + Cmi;
    dnew = GG'*dnorm;
end

% KEY LOOP OVER DIFFERENT MODELS (EACH OBTAINED FROM THE SAME EVENT KERNELS)
for ip = 1:nump
    stip = sprintf(stfm,ip);

    % OPTIONS FOR INVERTING HESSIAN:
    % (1) truncated singular value decomposition
    % (2) cross-validation
    if TSVD == 1
        pmax = pmax_vec(ip);
        disp(sprintf('%i out of %i : pmax = %i',ip,nump,pmax));
        
        %[mu_all,rss,f_r_ss] = tsvd(dnorm,H,pinds);   % KEY: TSVD
        mu = mu_all(:,pmax);        % KEY vector

        stplab = sprintf('pmax = %i',pmax);
        if 0 == 1
            figure; nr=2; nc=2;
            subplot(nr,nc,1);
            plot(pinds,s,'.',pinds(pmax),s(pmax),'pr','markersize',20,'linewidth',2);
            grid on; xlabel(xlab1); ylabel(ylab1); title({dir2lab,stplab});
            subplot(nr,nc,2);
            semilogy(pinds,s,'.',pinds(pmax),s(pmax),'pr','markersize',20,'linewidth',2);
            grid on; xlabel(xlab1); ylabel(ylab1); title(stplab);
            subplot(nr,nc,3);
            plot(pinds,rss,'.',pinds(pmax),rss(pmax),'pr','markersize',20,'linewidth',2);
            grid on; xlabel(xlab2); ylabel(ylab2); title(stplab);
            subplot(nr,nc,4);
            semilogy(pinds,rss,'.',pinds(pmax),rss(pmax),'pr','markersize',20,'linewidth',2);
            grid on; xlabel(xlab2); ylabel(ylab2); title(stplab);
            orient tall, wysiwyg
        end

    else
        if inew == 1
            midlampwr = log10(trace(Hnew))/2;
        else
            midlampwr = log10(trace(H))/2;
        end
        
        % regularization choices
        %minlampwr = -3; maxlampwr = 3;
        minlampwr = midlampwr - 2; maxlampwr = midlampwr + 2;
        numlam = 100;
        lampwr = linspace(minlampwr,maxlampwr,numlam);
        lamvec = 10.^lampwr;
        
        if inew == 1
            [f_p, rss, mss, Gvec, Fvec, dof, kap, iL, iGCV, iOCV] = ridge_carl(dnew,Hnew,lamvec);
        else
            [f_p, rss, mss, Gvec, Fvec, dof, kap, iL, iGCV, iOCV] = ridge_carl(dnorm,H,lamvec);
        end
        
        ipick = input(' Enter 0 for GCV, 1 for OCV, 2 for L-curve : ');
        switch ipick
            case 0, lam = lamvec(iGCV);
            case 1, lam = lamvec(iOCV);
            case 2, lam = lamvec(iL);
        end
        
        % KEY vector
        if inew == 1
            mu = inv(Hnew'*Hnew + lam^2*eye(3*nsrc,3*nsrc))*Hnew'*dnew;
        else
            mu = inv(H'*H + lam^2*eye(nsrc,nsrc))*H'*dnorm;
        end
        
        stplab = sprintf('lambda = %.3f',lam);
        
        %Hinv = inv(H);
        %mu = Hinv * dnorm;
    end
    
    % testing H + I
    %Hinv = inv(H);
    %mu = Hinv * dnorm;

    % KEY: solve for the STRUCTURE model update
    dm_str = zeros(nmod_str,1);
    if INV_STRUCT == 1
        if inew == 1
            dm_str = cov_beta' .* (Pstr * mu(1:nsrc));
        else
            dm_str = cov_beta' .* (transpose(Gstr) * mu);
        end
    end
    disp(' norm, min, max of structure update:');
    norm_dm_str = sum( dm_str .* dm_str ./ cov_beta' )
    min_dm_str = min( dm_str )
    max_dm_str = max( dm_str )
    
    % KEY: solve for the SOURCE model update (ts1, ts2, ..., xs1, xs2, ..., ys1, ys2, ...)
    dm_src = zeros(nmod_src,1);
    if INV_SOURCE == 1
        %mtemp = repmat(mu',3,1); mu_expand = mtemp(:);
        %wtemp = repmat(ws',3,1); ws_expand = wtemp(:);
        %dm_src = -cov_src' .* grad_src' .* ws_expand .* mu_expand;  % SIGN OR NOT?
        
        if inew == 1
            dm_src = cov_src' .* (Psrc * mu(nsrc+1:end));
        else
            dm_src = cov_src' .* (transpose(Gsrc) * mu);
        end
    end
    disp(' norm, min, max of source update:');
    norm_dm_src = sum( dm_src .* dm_src ./ cov_src' )
    min_dm_src = min( dm_src )
    max_dm_src = max( dm_src )
    
    %-------------------

    % source parameter updates
    if INV_SOURCE == 1

        %for ii = 1:nsrc
        %    i0 = (ii-1)*3;
        %    m_src_xs_new(ii) = m_src_xs(ii) + dm_src(i0+1);
        %    m_src_ys_new(ii) = m_src_ys(ii) + dm_src(i0+2);
        %    m_src_ts_new(ii) = m_src_ts(ii) + dm_src(i0+3);
        %end
        
        m_src_ts_new = m_src_ts + dm_src(inds_ts - nmod_str);
        m_src_xs_new = m_src_xs + dm_src(inds_xs - nmod_str);
        m_src_ys_new = m_src_ys + dm_src(inds_ys - nmod_str);
        
        m_src_ts_d_new = m_src_ts_new - m_src_ts_dat;
        m_src_xs_d_new = m_src_xs_new - m_src_xs_dat;
        m_src_ys_d_new = m_src_ys_new - m_src_ys_dat;

%         idisplay = 0;   % display information on source parameters
%         if idisplay == 1
%             disp('  '); disp('Source parameter updates:');
%             for ii = 1:nsrc
%                 i0 = (ii-1)*3;
%                 disp(sprintf('  Event %2i : (xs = %7.1f m,  ys = %7.1f m,  ts = %5.2f s )',...
%                    ii,dm_src(i0+1),dm_src(i0+2),dm_src(i0+3) ));
%             end
%             disp('  '); disp('Current parameter errors:');
%             for ii = 1:nsrc
%                 disp(sprintf('  Event %2i : (xs = %7.1f m,  ys = %7.1f m,  ts = %5.2f s )',...
%                    ii, m_src_xs_d(ii), m_src_ys_d(ii),m_src_ts_d(ii)));
%             end
%             disp('  '); disp('New source parameter errors:');
%             for ii = 1:nsrc
%                 disp(sprintf('  Event %2i : (xs = %7.1f m,  ys = %7.1f m,  ts = %5.2f s )',...
%                    ii, m_src_xs_d_new(ii), m_src_ys_d_new(ii),m_src_ts_d_new(ii)));
%             end
%         end
    else
        m_src_ts_d_new = m_src_ts_d;
        m_src_xs_d_new = m_src_xs_d;
        m_src_ys_d_new = m_src_ys_d;
        
        m_src_ts_new = m_src_ts;
        m_src_xs_new = m_src_xs;
        m_src_ys_new = m_src_ys;
    end

    % compute norms of errors in new source parameters
    %norm_ts = sum(m_src_ts_d_new.^2 ./ cov_src(3:3:nsrc*3)');
    %norm_xs = sum(m_src_xs_d_new.^2 ./ cov_src(1:3:nsrc*3)');
    %norm_ys = sum(m_src_ys_d_new.^2 ./ cov_src(2:3:nsrc*3)');
    norm_ts = sum(m_src_ts_d_new.^2 ./ cov_model(inds_ts));
    norm_xs = sum(m_src_xs_d_new.^2 ./ cov_model(inds_xs));
    norm_ys = sum(m_src_ys_d_new.^2 ./ cov_model(inds_ys));
    source_error_norm_mat(ip,[2:5]) = [norm_xs norm_ys norm_ts sum([norm_xs norm_ys norm_ts])];
    
    % structure parameter updates
    if INV_STRUCT == 1
        m_str_B_new = m_str_B + dm_str;   % updated B
    else
        m_str_B_new = m_str_B;
    end

    % convert to beta, c, kappa, mu, rho
    m_str_beta_new  = beta0 * exp( m_str_B_new );
    m_str_rho_new   = m_str_rho;
    m_str_kappa_new = m_str_kappa;
    m_str_mu_new    = m_str_rho_new .* m_str_beta_new.^2;

    %disp(' norm of new structure model:');
    %sum( m_str_B_new .* m_str_B_new ./ cov_beta' )
    
    % PLOT structure model update and new model
    if and(INV_STRUCT == 1, ifigure == 1)
        tlabs = {{dir2lab,'data'},'m00',['dm (' stplab ')'],['m01 (' stplab ')']};
        cmax = 0.1;
        figure; nr=3; nc=2;
        for kk = 1:4
            switch kk
                case 1, mplot = m_dat_str_B;
                case 2, mplot = m_str_B;
                case 3, mplot = dm_str;
                case 4, mplot = m_str_B_new;
            end
            [X,Y,Z] = griddataXB(lon,lat,mplot,100,'nearest');
            mnorm = sum( mplot.^2 ./ cov_beta' );

            subplot(nr,nc,kk); hold on;
            pcolor(X,Y,Z); shading interp;
            caxis([-1 1]*cmax); colormap(seis); axis equal, axis tight
            plot(rlon,rlat,'k.','markersize',16);
            plot(slon,slat,'p','markersize',18,'markeredgecolor','k','linewidth',2);
            title([tlabs{kk} sprintf(', norm = %.2f',mnorm)]);
        end
        
        % TSVD plots
        if TSVD == 1
            subplot(nr,nc,5);
            semilogy(pinds,s,'.',pinds(pmax),s(pmax),'pr','markersize',20,'linewidth',2);
            grid on; xlabel(xlab1); ylabel(ylab1); title(stplab);
            subplot(nr,nc,6);
            semilogy(pinds,rss,'.',pinds(pmax),rss(pmax),'pr','markersize',20,'linewidth',2);
            grid on; xlabel(xlab2); ylabel(ylab2); title(stplab);
        end
        
        orient tall, wysiwyg
    end

    % write EVERY pmax updated model to file
    % NOTE: if you only want one model, then it will have index 001
    if iwrite == 1
        
        % KEY: the directories differ, depending on whether you are writing
        % out a single model or a set of models
        if COMPUTE_KERNELS
            dir3 = dir2;
        else
            % make directory
            pdir = ['run_p' stip];
            dir3 = [dir2 pdir '/'];
            mkdir(dir3)
        end
        
        % save pertinent Matlab variables
        save([dir3 'wave2d_subspace_matlab_m' stm],...
            'dnorm','H','Hsrc','Hstr',...
            'mu','INV_SOURCE','INV_STRUCT');

        % write structure model
        wave2d_write_str([dir3 'structure_syn_m' stm '.dat'],m_str_lon,m_str_lat,...
            m_str_kappa_new,m_str_mu_new,m_str_rho_new,m_str_B_new);
        
        % write source model
        wave2d_write_src([dir3 'src_syn_m' stm '.dat'],m_src_lon,m_src_lat,...
            m_src_ts_new,m_src_xs_new,m_src_ys_new,m_src_ts_d_new,m_src_xs_d_new,m_src_ys_d_new);

        % write out pmax
        if TSVD == 1
            fid = fopen([dir3 'pmax_m' stm '.dat'],'w');
            fprintf(fid,'%i\n',pmax);   
            fclose(fid);
        else
            fid = fopen([dir3 'lambda_m' stm '.dat'],'w');
            fprintf(fid,'%16.8e\n',lam);   
            fclose(fid);
        end
        
        % write out the balance of the gradient norms
        fid = fopen([dir3 'norm_gradient_m' stm '.dat'],'w');
        fprintf(fid,'Balance of the gradients for each event\n');
        fprintf(fid,'%16s%16s%16s%16s\n','STR','SRC','TOTAL','SRC/STR');
        for ii = 1:nsrc
            fprintf(fid,'%16.6e%16.6e%16.6e%16.4f\n',...
                Hstr(ii,ii), Hsrc(ii,ii), H(ii,ii), Hsrc(ii,ii)/Hstr(ii,ii) );
        end
        fprintf(fid,'%48s%16.4f\n','MEAN -->',mean( diag(Hsrc)./diag(Hstr)) );
        fclose(fid);
        
    end
    
end  % LOOP OVER ip

% norm in intial source errors
%norm_xs = sum(m_src_xs_d.^2 ./ cov_src(1:3:nsrc*3)');
%norm_ys = sum(m_src_ys_d.^2 ./ cov_src(2:3:nsrc*3)');
%norm_ts = sum(m_src_ts_d.^2 ./ cov_src(3:3:nsrc*3)');
norm_ts = sum(m_src_ts_d.^2 ./ cov_model(inds_ts));
norm_xs = sum(m_src_xs_d.^2 ./ cov_model(inds_xs));
norm_ys = sum(m_src_ys_d.^2 ./ cov_model(inds_ys));
source_error_norm_old = [0 norm_xs norm_ys norm_ts sum([norm_xs norm_ys norm_ts])];

disp('   ');
disp(' NORMS of the errors in the old source vectors (xs, ys, ts) :');
disp(source_error_norm_old);
disp(' NORMS of the errors in the new source vectors (xs, ys, ts) :');
disp('    pmax         xs        ys        ts      total');
disp(source_error_norm_mat);

% save pertinent Matlab variables
if iwrite == 1
    if COMPUTE_KERNELS
        ofile = [dir2 'wave2d_subspace_matlab_pmax'];
    else
        ofile = [dir2 'wave2d_subspace_matlab_all'];
    end
    
    if TSVD == 0
        save(ofile,'dnorm','lamvec','lam',...
            'f_p','rss','mss','Gvec','Fvec','dof','kap','iL','iGCV','iOCV',...
            'source_error_norm_mat','source_error_norm_old',...
            'INV_SOURCE','INV_STRUCT');
    else   
        save(ofile,'dnorm','U','S','V','pinds','mu_all','rss','f_r_ss',...
            'source_error_norm_mat','source_error_norm_old',...
            'INV_SOURCE','INV_STRUCT');
    end
end

%======================================================================
% EXTRA CODE FOR TESTING

%------------------------------------------------
% extra code for plotting misfit values

if 0==1
   clear
   
   run0 = 6550; close all
   for h = 1:4

       strun0 = num2str(run0);
       stm = sprintf(stfm,h);
       dir_run = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_OUTPUT/';
       idir = [dir_run 'run_' strun0 '/READ_IN/model_m' stm '/'];
       load([idir 'wave2d_subspace_matlab_all']);

       clear data_norm model_norm chi_total pmax_vec
       for ii = 1:25
          stp = sprintf(stfm,ii);
          dir1 = [ idir 'run_p' stp '/'];
          chi_file = [dir1 'data_norm.dat'];
          if exist(chi_file)
              data_norm(ii) = load(chi_file);
              model_norm(ii) = load([dir1 'model_norm.dat']);
              chi_total(ii) = load([dir1 'chi.dat']);
              %chi_total(ii) = load([dir1 'summed_chi_all.dat']);
              pmax_vec(ii) = load([dir1 'pmax_m' stm '.dat']);
          end
       end

       %disp('  '); disp('     pmax  2*S(m(pmax)) S-data   S-model');
       %disp([pinds 2*chi_total' data_norm' model_norm']);

       % semilogy plots or not?
       figure; nr=2; nc=2;
       xlab1 = 'p, singular value index';
       xlab2 = 'p, singular value truncation index';
       ylab1 = 'singular value';
       ylab2 = 'misfit : dot[ d - G*dm(p), d - G*dm(p) ]';

       subplot(nr,nc,1); semilogy(pinds, diag(S),'.','markersize',22); grid on;
       xlabel(xlab1); ylabel(ylab1);
       subplot(nr,nc,3); semilogy(pinds, rss,'.','markersize',22); grid on;
       xlabel(xlab2); ylabel(ylab2);

       subplot(nr,nc,2);
       plot( pmax_vec,2*chi_total,'r.',pmax_vec,data_norm,'b.','markersize',22);
       grid on; legend('2 S(m(pmax))','chi-data-norm');
       xlabel(xlab2); title(['RUN ' strun0 ' - H' stm]);
       subplot(nr,nc,4);
       plot( pmax_vec,model_norm,'g.','markersize',22);
       grid on; legend('chi-model-norm','location','northeast');
       xlabel(xlab2); title(['RUN ' strun0 ' - H' stm]);

       orient tall, wysiwyg
   end
   
end

%----------------------------------------------

if 0==1
    temp = load('/net/denali/scratch1/carltape/OUTPUT/run_7000/structure_dat.dat');
    %temp = load('/net/denali/scratch1/carltape/OUTPUT/run_7000/READ_IN/structure_dat.dat');
    lon = temp(:,1); lat = temp(:,2); beta = temp(:,7);

    [X,Y,Z] = griddataXB(lon,lat,beta,100,'nearest');
    figure; cmax = 0.1; hold on;
    pcolor(X,Y,Z); shading interp;
    caxis([-1 1]*cmax); colormap(seis);
end

%----------------------------------------------

% check the relative dimensions of matrices and plot a schematic figure
if 0==1
    S = 5;              % number of sources
    Cstr = 1;           % number of structure parameter classes
    Csrc = 2;           % number of source parameter classes
    C = Cstr + Csrc;
    Mstr1 = 20;         % number of structure parameters
    Msrc1 = 1;          % number of source parameters in class 1 (origin time)
    Msrc2 = 2;          % number of source parameters in class 2 (xy, ys)
    Mstr = Mstr1;
    Msrc = S*(Msrc1 + Msrc2);
    M = Mstr + Msrc;    % number of rows in P
    N = C*S;            % number of columns in P
    
    P = zeros(M,N);
    P(1:Mstr1,1:S) = rand(Mstr1,S);
    for s = 1:S
       P(Mstr1+s,S+s) = rand; 
       P(Mstr1+S+s,2*S+s) = rand; 
       P(Mstr1+2*S+s,2*S+s) = rand; 
    end
    
    % fill G
    G = zeros(S,M);
    Gstr = rand(S,Mstr);
    Gsrc = zeros(S,Msrc);
    if 1==1
        Gsrc = repmat(diag(rand(S,1)),1,Msrc1+Msrc2);
    else
        itemp = (Msrc1 + Msrc2)*[1:S]';
        indmat2 = [ [1 ; 1+itemp(1:end-1)] itemp ];
        for isrc = 1:S
            inds = indmat2(isrc,1) : indmat2(isrc,2);
            %Gsrc(isrc,inds) = -ws(isrc) * grad_src(inds);
            Gsrc(isrc,inds) = rand(1,length(inds));
        end
    end
    G = [Gstr Gsrc];
   
    Cm  = diag(rand(M,1));
    Cmi = P'*Cm*P;
    GG  = G*Cm*P;
    
    d = zeros(S,1);
    Hnew = GG'*GG + Cmi;
    dnew = GG'*d;
    mu = Hnew*dnew;
    whos d Hnew dnew mu
    
    figure; nr=3; nc=4; msize=3;
    %subplot(nr,nc,2); spy(G,msize); title(sprintf('G is %i by %i',size(G)));
    %subplot(nr,nc,3); spy(G',msize); title(sprintf('G^T is %i by %i',size(G')));
    %subplot(nr,nc,4); spy(G*G',msize); title(sprintf('G G^T is %i by %i',size(G*G')));
    subplot(nr,nc,1); spy(G,msize); title(sprintf('G is %i by %i',size(G)));
    subplot(nr,nc,2); spy(Cm,msize); title(sprintf('Cm is %i by %i',size(Cm)));
    subplot(nr,nc,3); spy(G',msize); title(sprintf('G^T is %i by %i',size(G')));
    subplot(nr,nc,4); spy(G*Cm*G',msize); title(sprintf('G Cm G^T is %i by %i',size(G*Cm*G')));
    
    subplot(nr,nc,5); spy(P',msize); title(sprintf('P^T is %i by %i',size(P')));
    subplot(nr,nc,6); spy(Cm,msize); title(sprintf('Cm is %i by %i',size(Cm)));
    subplot(nr,nc,7); spy(P,msize); title(sprintf('P is %i by %i',size(P)));
    subplot(nr,nc,8); spy(Cmi,msize); title(sprintf('P^T Cm P is %i by %i',size(Cmi)));
    
    subplot(nr,nc,9); spy(G,msize); title(sprintf('G is %i by %i',size(G)));
    subplot(nr,nc,10); spy(Cm,msize); title(sprintf('Cm is %i by %i',size(Cm)));
    subplot(nr,nc,11); spy(P,msize); title(sprintf('P is %i by %i',size(P)));
    subplot(nr,nc,12); spy(GG,msize); title(sprintf('G Cm P is %i by %i',size(GG)));
    orient tall, wysiwyg
    
    % projection for data
    S = 5;
    N = 15;
    Ns = N/S;
    
    P = zeros(S,N);
    for s = 1:S
       P(s, Ns*(s-1) + [1:Ns]) = rand; 
    end
    d = randn(N,1);
    C = rand(N,N);
    Cblock = P'*P;
    Cdiag = diag(rand(N,1));
    
    figure; nr=4; nc=4;
    subplot(nr,nc,1); spy(P); title(sprintf('P is %i by %i',size(P))); 
    subplot(nr,nc,2); spy(P'); title(sprintf('P^T is %i by %i',size(P')));
    subplot(nr,nc,3); spy(P*P'); title(sprintf('P*P^T is %i by %i',size(P*P')));
    
    subplot(nr,nc,5); spy(P'); title(sprintf('P^T is %i by %i',size(P')));
    subplot(nr,nc,6); spy(P); title(sprintf('P is %i by %i',size(P)));
    subplot(nr,nc,7); spy(P'*P); title(sprintf('P^T*P is %i by %i',size(P'*P)));
    
    subplot(nr,nc,9); spy(C'); title(sprintf('C^T is %i by %i',size(C')));
    subplot(nr,nc,10); spy(P'*P); title(sprintf('P^T*P is %i by %i',size(P'*P)));
    subplot(nr,nc,11); spy(C); title(sprintf('C is %i by %i',size(C)));
    subplot(nr,nc,12); spy(C'*P'*P*C); title(sprintf('C^T*P^T*P*C is %i by %i',size(C'*P'*P*C)));
  
    subplot(nr,nc,13); spy(Cdiag'); title(sprintf('C^T is %i by %i',size(Cdiag')));
    subplot(nr,nc,14); spy(P'*P); title(sprintf('P^T*P is %i by %i',size(P'*P)));
    subplot(nr,nc,15); spy(Cdiag); title(sprintf('C is %i by %i',size(Cdiag)));
    subplot(nr,nc,16); spy(Cdiag'*P'*P*Cdiag); title(sprintf('C^T*P^T*P*C is %i by %i',size(Cdiag'*P'*P*Cdiag)));
    orient tall, wysiwyg
    %print(gcf,'-dpsc', 'test')

    figure; nr=4; nc=4;
    subplot(nr,nc,1); spy(P); title(sprintf('P is %i by %i',size(P))); 
    subplot(nr,nc,2); spy(C); title(sprintf('C is %i by %i',size(C)));
    subplot(nr,nc,3); spy(d); title(sprintf('d is %i by %i',size(d)));
    subplot(nr,nc,4); spy(P*C*d); title(sprintf('P*C*d is %i by %i',size(P*C*d)));
    
    subplot(nr,nc,5); spy(P); title(sprintf('P is %i by %i',size(P))); 
    subplot(nr,nc,6); spy(Cblock); title(sprintf('C is %i by %i',size(Cblock)));
    subplot(nr,nc,7); spy(d); title(sprintf('d is %i by %i',size(d)));
    subplot(nr,nc,8); spy(P*Cblock*d); title(sprintf('P*C*d is %i by %i',size(P*Cblock*d)));
    
    subplot(nr,nc,9); spy(P); title(sprintf('P is %i by %i',size(P))); 
    subplot(nr,nc,10); spy(Cdiag); title(sprintf('C is %i by %i',size(Cdiag)));
    subplot(nr,nc,11); spy(d); title(sprintf('d is %i by %i',size(d)));
    subplot(nr,nc,12); spy(P*Cdiag*d); title(sprintf('P*C*d is %i by %i',size(P*Cdiag*d)));
    
    subplot(nr,nc,14); spy(P*Cdiag); title(sprintf('P*C is %i by %i',size(P*Cdiag))); 
    subplot(nr,nc,15); spy(Cdiag*d); title(sprintf('C*d is %i by %i',size(Cdiag*d)));
    orient tall, wysiwyg
    %print(gcf,'-dpsc', 'test')
    
    break
end

%---------------------------------------------

if 0==1
    nk = 100;
    kvec = 10.^linspace(-4,4,nk);
    
    figure; nr=2; nc=1;
    for jj=1:2
        switch jj
            case 1, ndata = 5; nparm = 10;
            case 2, ndata = 10; nparm = 5;
        end
        Cm = diag(rand(nparm,1));
        G = rand(ndata,nparm);
        d = rand(ndata,1);
        
        norm_dm = zeros(nk,1);
        for kk = 1:nk
            C = Cm * kvec(kk);
            norm_dm(kk) = norm( C * G' * inv(G*C*G') * d );
        end
        subplot(nr,nc,jj);
        semilogx(kvec, norm_dm, '.');
    end
end

%======================================================================
