%
% wave2d_sigmas.m
% CARL TAPE, 20-Jan-2010
% printed xxx
% 
% Produces Gaussian errors using the Matlab function randn.
%   randn(1,N) generates an N x 1 column vector, the values having a
%   standard deviation of 1.
% 
% calls plot_histo.m
% called by xxx
%

format short
format compact
close all
clear
clc

igauss_diff = 0;
isigma_source = 1;
isigma_struct = 0;
isigma_data = 0;

iwrite = 0;

%odir = '/net/denali/home2/carltape/wave2d/2d_adjoint/INPUT/';
odir = '/home/carltape/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work/INPUT/';

% labels for source parameters
xlabs = {'origin time (s)','x position (m)','y position (m)'};
xlabs2 = {'ts (s)','xs (m)','ys (m)'};

NLOCAL = 40000;
NPARM_SOURCE = 3;

%---------------------------------------------------------
% demonstrate the analytical expression for a linear combination of two Gaussians

if igauss_diff == 1
    
    % input parameters
    N = 1000;
    sig1 = 2; mu1 = 4;
    sig2 = 6; mu2 = -20;
    A = 0.5; B = -3;
    
    g1 = mu1 + sig1 * randn(N,1);
    g2 = mu2 + sig2 * randn(N,1);
    dmat = zeros(N,N);
    for ii=1:N
        for jj=1:N
            % linear sum of every pair of samples
            dmat(ii,jj) = A*g1(ii) + B*g2(jj);
        end
    end
    g1g2 = dmat(:);

    % KEY: analytical expressions for combined mu and sig
    sigs = [sig1 sig2 sqrt(A^2*sig1^2 + B^2*sig2^2)];
    mus  = [mu1 mu2 A*mu1+B*mu2];

    figure; nr=3; nc=1; ylims = [0 0.15];
    subplot(nr,nc,1); sig = sigs(1); mu = mus(1);
    Nbin = plot_histo(g1,mu + [-4*sig : sig/4 : 4*sig]);
    grid on; xlabel(sprintf('g1: mean %.2f sigma %.2f',mu,sig)); ylim(ylims);

    subplot(nr,nc,2); sig = sigs(2); mu = mus(2);
    Nbin = plot_histo(g2,mu + [-4*sig : sig/4 : 4*sig]);
    grid on; xlabel(sprintf('g2: mean %.2f sigma %.2f',mu,sig)); ylim(ylims);

    subplot(nr,nc,3); sig = sigs(3); mu = mus(3);
    Nbin = plot_histo(g1g2,mu + [-4*sig : sig/4 : 4*sig]);
    grid on;  xlabel(sprintf('A*g1 + B*g2: mean %.2f sigma %.2f (A = %.2f, B = %.2f)',mu,sig,A,B)); ylim(ylims);
    orient tall, wysiwyg;
    
end

%---------------------------------------------------------

if isigma_source == 1
    ifig = 0;
    
    %N = 10000; RUNS = 1;
    N = 25; RUNS = 1000;
    
    % KEY COMMANDS
    %sigs = [0.5 2000 2000];    % sigmas to approximate the GJI2007 perturbations
    sigs = [0.2 1400 1400];
    
    % the STD for the difference (or sum) of two samples
    sigssum = sqrt( 2*sigs.^2 )
    
    cov_model = (N*NPARM_SOURCE) * repmat(sigs.^2,N,1);
    
    % generate RUNS samples
    ts_samples = sigs(1)*randn(N,RUNS);
    xs_samples = sigs(2)*randn(N,RUNS);
    ys_samples = sigs(3)*randn(N,RUNS);
    
    % compute norms
    norm2_ts = sum( ts_samples.^2 / (N*NPARM_SOURCE*sigs(1)^2), 1)';
    norm2_xs = sum( xs_samples.^2 / (N*NPARM_SOURCE*sigs(2)^2), 1)';
    norm2_ys = sum( ys_samples.^2 / (N*NPARM_SOURCE*sigs(3)^2), 1)';
    norm2_all = [norm2_ts norm2_xs norm2_ys];
    norm2_tot = sum(norm2_all,2);
    
    figure; nr=2; nc=2;
    ncen = 1/3;
    for ii=1:3
        subplot(nr,nc,ii);
        Nbin = plot_histo(norm2_all(:,ii),[0:ncen/8:2*ncen]);
        grid on; xlabel(['Norm of ' xlabs{ii} ' for each sample']);
        ylim([0 0.3]);
    end
    ncen = 1; subplot(nr,nc,4);
    Nbin = plot_histo(norm2_tot,[0:ncen/8:2*ncen]);
    grid on; xlabel('Norm of each sample');
    ylim([0 0.3]);
    orient tall, wysiwyg;
    
    % pick a subset of samples whose overall norms are close to 1.0 and
    % whose and norm-parts are close to 0.33
    i1 = find( abs( log(norm2_tot / 1.0 )) < 0.1 );
    i2 = find( abs( log(norm2_ts / (1/3) )) < 0.1 );
    i3 = find( abs( log(norm2_xs / (1/3) )) < 0.1 );
    i4 = find( abs( log(norm2_ys / (1/3) )) < 0.1 );
    iset = mintersect(i1,i2,i3,i4);
    nset = length(iset)
    
    if ifig==1
        for kk=1:nset
            jj = iset(kk);
            Ymat = [ts_samples(:,jj) xs_samples(:,jj) ys_samples(:,jj)];
            Nvec = [norm2_all(jj,:) norm2_tot(jj)];

            figure; nr=3; nc=1;
            for ii = 1:3 
                sig = sigs(ii);
                yvec = Ymat(:,ii);
                subplot(nr,nc,ii);
                Nbin = plot_histo(yvec,[-5*sig : sig/4 : 5*sig]);
                grid on; xlabel(['Perturbation in ' xlabs{ii}]); ylim([0 0.2]);
                title(sprintf('%s -- sigma = %.3f -- Norm of %i parameters = %.4f',...
                    xlabs{ii},sig,N,Nvec(ii)));
            end
            orient tall, wysiwyg;
        end
    end
    
    % from the set, pick two samples whose DIFFERENCE is not an outlier
    % --> the DIFFERENCE represents the net perturbation needed
    n = 0;
    Imat = zeros(n,3);
    Nmat = zeros(n,3);
    for kk=1:nset
        ik = iset(kk);
        Kmat = [ts_samples(:,ik) xs_samples(:,ik) ys_samples(:,ik)];
        for jj=kk+1:nset
            ij = iset(jj);
            Jmat = [ts_samples(:,ij) xs_samples(:,ij) ys_samples(:,ij)];
            Rmat = Kmat - Jmat;
            
            n = n+1;
            % NOTE siggssum instead of sigs
            Nmat(n,1) = sum( Rmat(:,1).^2 / (N*NPARM_SOURCE*sigssum(1)^2) );
            Nmat(n,2) = sum( Rmat(:,2).^2 / (N*NPARM_SOURCE*sigssum(2)^2) );
            Nmat(n,3) = sum( Rmat(:,3).^2 / (N*NPARM_SOURCE*sigssum(3)^2) );
            Imat(n,:) = [n ik ij];
        end
    end
    Nmat_tot = sum(Nmat, 2);
    %figure; plot(Nmat_tot,'.')
    
    % pick a pair whose norm of the DIFFERENCE (and its parts) is not an outlier
    i1 = find( abs( log(Nmat_tot / 1.0 )) < 0.1 );
    i2 = find( abs( log(Nmat(:,1) / (1/3) )) < 0.1 );
    i3 = find( abs( log(Nmat(:,2) / (1/3) )) < 0.1 );
    i4 = find( abs( log(Nmat(:,3) / (1/3) )) < 0.1 );
    ipairs = mintersect(i1,i2,i3,i4);
    npair = length(ipairs);
    %Imat(ipairs,:)
    
    % plot the pairs
    for ii=1:npair
        inds = Imat(ipairs(ii),:);
        i1 = inds(2);
        i2 = inds(3);
        
        g1 = [ts_samples(:,i1) xs_samples(:,i1) ys_samples(:,i1)];
        g2 = [ts_samples(:,i2) xs_samples(:,i2) ys_samples(:,i2)];
        g1g2 = g2 - g1;
        
        % event 5
        d5 = sqrt( g1g2(5,2)^2 + g1g2(5,3)^2 );
        t5 = g1g2(5,1);
        disp(sprintf('%6i%6i%6.2f s%10.1f m',i1,i2,t5,d5));
        
        if ifig == 1
            figure; nr=3; nc=3;
            for kk=1:3
                subplot(nr,nc,kk); sig = sigs(kk);
                Nbin = plot_histo(g1(:,kk),[-4*sig : sig/4 : 4*sig]);
                title(sprintf('sample %i (%.2f)',i1,max(abs(g1(:,kk)))/sig)); xlabel(xlabs2{kk});
                subplot(nr,nc,3+kk); sig = sigs(kk);
                Nbin = plot_histo(g2(:,kk),[-4*sig : sig/4 : 4*sig]);
                title(sprintf('sample %i (%.2f)',i2,max(abs(g2(:,kk)))/sig)); xlabel(xlabs2{kk});

                subplot(nr,nc,6+kk); sig = sigssum(kk);
                Nbin = plot_histo(g1g2(:,kk),[-4*sig : sig/4 : 4*sig]);
                title(sprintf('sample %i - %i (%.2f)',i2,i1,max(abs(g1g2(:,kk)))/sig)); xlabel(xlabs2{kk});

            end
            fontsize(8), orient landscape, wysiwyg;
        end
    end
    
%     nr=3; nc=1;
%     for kk = 1:RUNS
%         kk
% 
%         if ifig==1, figure; end
%         Ymat = zeros(N,3);
%         for ii = 1:3            % loop over three parameters
%             sig = sigs(ii);
%             yvec = sig*randn(N,1);
%             Ymat(:,ii) = yvec;
%             
%             if ifig==1
%                 subplot(nr,nc,ii);
%                 Nbin = plot_histo(yvec,[-5*sig : sig/4 : 5*sig]);
%                 grid on; xlabel(['Perturbation in ' xlabs{ii}]);
%                 ylim([0 0.2]);
%             end
%         end
%     end
%     
%     for kk = 1:NRUNS
%         
%         % pick the sample based on
%         %   (1) the ith event perturbation
%         %   (2) the overall norm (close to 1)
%         %   (3) the STD of the norms of each parameter (low)
%         
%         % ith event
%         %irec = 5; dloc_min = 4000; dotime_min = 0.8;
%         irec = 5; dloc_min = 4000; dotime_min = 0.2;
%         misloc = sqrt( Ymat(irec,2)^2 + Ymat(irec,3)^2 );
%         otime = Ymat(irec,1);
%         
%         % overall norm
%         norms = sum( Ymat .* Ymat ./ cov_model);
%         
%         if and( and( misloc > dloc_min, abs(otime) > dotime_min ), ...
%                 and( abs( sum(norms) - 1 ) <= 0.05, std(norms) <= 0.02) )
%             disp(' Event 5:');
%             disp([' Mislocation (m) : ' num2str(misloc) ]);
%             disp([' Origin time perturbation (s) : ' num2str(otime) ]);
%             
%             figure;
%             for ii = 1:3 
%                 sig = sigs(ii); yvec =Ymat(:,ii);
%                 subplot(nr,nc,ii);
%                 Nbin = plot_histo(yvec,[-5*sig : sig/4 : 5*sig]);
%                 grid on; xlabel(['Perturbation in ' xlabs{ii}]); ylim([0 0.2]);
%                 title(sprintf('%s -- sigma = %.3f -- Norm of %i parameters = %.4f',...
%                     xlabs{ii},sig,N,norms(ii)));
%             end
%             orient tall, wysiwyg;
%             break
%         end
%         if ifig==1, orient tall, wysiwyg; end
%     end
    
    if iwrite==1
        disp('Are you sure you want to write two samples to files?');
        i1 = input('  Enter index of initial model: '); 
        i2 = input('  Enter index of target model: '); 
        Y1 = [ts_samples(:,i1) xs_samples(:,i1) ys_samples(:,i1)];
        Y2 = [ts_samples(:,i2) xs_samples(:,i2) ys_samples(:,i2)];
        YD = Y2 - Y1;
        
        % check norms
        norm1 = [  sum( Y1(:,1).^2 / (N*NPARM_SOURCE*sigs(1)^2) ) 
                   sum( Y1(:,2).^2 / (N*NPARM_SOURCE*sigs(2)^2) ) 
                   sum( Y1(:,3).^2 / (N*NPARM_SOURCE*sigs(3)^2) )  ]
        sum(norm1)
        norm2 = [  sum( Y2(:,1).^2 / (N*NPARM_SOURCE*sigs(1)^2) ) 
                   sum( Y2(:,2).^2 / (N*NPARM_SOURCE*sigs(2)^2) ) 
                   sum( Y2(:,3).^2 / (N*NPARM_SOURCE*sigs(3)^2) )  ]
        sum(norm1)
        normD = [  sum( YD(:,1).^2 / (N*NPARM_SOURCE*sigssum(1)^2) ) 
                   sum( YD(:,2).^2 / (N*NPARM_SOURCE*sigssum(2)^2) ) 
                   sum( YD(:,3).^2 / (N*NPARM_SOURCE*sigssum(3)^2) )  ]
        sum(normD)
        
        ofile = [odir 'events_txy_pert_initial.dat'];
        fid = fopen(ofile,'w');
        for ii = 1:N
            fprintf(fid,'%20.12e%20.12e%20.12e\n', Y1(ii,:));   
        end
        fclose(fid);
        
        ofile = [odir 'events_txy_pert_target.dat'];
        fid = fopen(ofile,'w');
        for ii = 1:N
            fprintf(fid,'%20.12e%20.12e%20.12e\n', Y2(ii,:));   
        end
        fclose(fid);
        
        ofile = [odir 'events_txy_pert_sigmas.dat'];
        fid = fopen(ofile,'w');
        fprintf(fid,'%20.12e%20.12e%20.12e\n',sigs);   
        fclose(fid);
    end
    
    %xmin = ans(1); xmax = ans(2);
    %xpt = linspace(xmin,xmax,100);
    %ypt = 1/sqrt(2*sig^2) * exp(-xpt.^2 / (2*sig^2) )
end

%---------------------------------------------------------
% see what the sine-like checkerboard maps look like as a distribution
% --> What is the sigma value that characterizes the distribution?

if isigma_struct == 1
    
    idir = '/net/denali/scratch1/carltape/OUTPUT_2/run_9100/';  % 7100
    temp = load([idir 'structure_dat.dat']); beta = temp(:,7);
    N = length(beta);
    
    sig0 = 0.1;
    sfac = [1 2 3 4 5];
    
    for kk = 1:length(sfac);
        sig = sig0 / sfac(kk);
    
        figure; nr=2; nc=1;
        subplot(nr,nc,1);
        Nbin = plot_histo(beta,[-5*sig : sig/4 : 5*sig]);
        grid on; xlabel(['Perturbation in ' xlabs{1}]);
        title('checkerboard with +/- 0.10');
        ylim([0 0.2]); %xlim([-0.5 0.5]);

        subplot(nr,nc,2);
        Nbin = plot_histo(sig*randn(N,1),[-5*sig : sig/4 : 5*sig]);
        grid on; xlabel(['Perturbation in ' xlabs{1}]);
        title(['sigma = ' sprintf('%.3f',sig)]);
        ylim([0 0.2]); %xlim([-0.5 0.5]);
        orient tall, wysiwyg
    end
    
    %---------------------
    
    % load the jacobian for constructing the "event gradient" from the event kernel
    lmesh_all = load([idir 'local_mesh.dat']);
    Ai = lmesh_all(:,9);    % local dA area element
    Atot = sum(Ai);         % total area
    
    % load the model covariance
    cov_model = load([idir 'cov_imetric_diagonal.dat']);
    cov_str   = cov_model(1:NLOCAL);
    
    sig0 = 0.05
    beta_const = sig0 * ones(N,1);
    beta_gaus  = sig0 * randn(N,1);
    
    % compute norms
    Mnorm0_gaus  = sum( beta_gaus.^2 ./ cov_str )
    Mnorm0       = sum( beta.^2 ./ cov_str )   
    Mnorm0_const = sum( beta_const.^2 ./ cov_str )
    
    %---------------------
    
    figure; nr=3; nc=1;
    bins = [-5*sig0 : sig0/4 : 5*sig0];
    
    subplot(nr,nc,1);
    Nbin = plot_histo(sig0*randn(N,1),bins);
    grid on; xlabel(['Perturbation in ' xlabs{1}]);
    title(sprintf('Gaussian with sigma = %.3f; norm (m^T C^{-1} m) = %.3f',...
        sig0,Mnorm0_gaus)); ylim([0 0.2]);
    
    subplot(nr,nc,2);
    Nbin = plot_histo(beta,bins);
    grid on; xlabel(['Perturbation in ' xlabs{1}]);
    title(sprintf('Checkerboard with amplitude 0.100; norm (m^T C^{-1} m) = %.3f',...
        sig0,Mnorm0)); ylim([0 0.2]);

    subplot(nr,nc,3);
    Nbin = plot_histo(beta_const,bins);
    grid on; xlabel(['Perturbation in ' xlabs{1}]);
    title(sprintf('Constant function with values = %.3f; norm (m^T C^{-1} m) = %.3f',...
        sig0,Mnorm0_const)); ylim([0 0.2]);
        
    orient tall, wysiwyg
    
    %---------------------------
    
    sigvec = 10.^linspace(-2,0,121);   % test sigma values
    for kk = 1:length(sigvec);
        sig = sigvec(kk);
        Mnorm(kk)       = (1/Atot) * sum(beta.^2 / sig^2 .* Ai);
        Mnorm_const(kk) = (1/Atot) * sum(beta_const.^2 / sig^2 .* Ai);
    end
    
    figure; nr=2; nc=2;
    st1 = ['checkerboard (PERT = 0.100)'];
    st2 = ['uniform (PERT = ' sprintf('%.3f',sig0) ')'];
    
    subplot(nr,nc,1);
    plot(sigvec, Mnorm,'b.', sigvec, Mnorm_const,'r.'); grid on;
    legend(st1,st2); xlabel(' sigma value'); ylabel(' Norm');
    
    subplot(nr,nc,2);
    loglog(sigvec, Mnorm,'b.', sigvec, Mnorm_const,'r.'); grid on;
    legend(st1,st2); xlabel(' sigma value'); ylabel(' Norm');

    subplot(nr,nc,3);
    plot(sigvec, Mnorm,'b.', sigvec, Mnorm_const,'r.'); grid on;
    legend(st1,st2); xlabel(' sigma value'); ylabel(' Norm');
    axis([0.03 0.12 0 2]);
    
    orient tall, wysiwyg
end

%---------------------------------------------------------
% see what the sine-like checkerboard maps look like as a distribution
% --> What is the sigma value that characterizes the distribution?

if isigma_data == 1
    
    %-------------------------------
    % load traveltime measurments from 2D code
    
    idir = '/net/denali/scratch1/carltape/OUTPUT_2/run_7000/';
    meas_all = load([idir 'measure_vec.dat']);
    dT_all = meas_all(:,2);
    Nmeas = length(dT_all);
    
    sig = 1.7; 
    figure; nr=2; nc=1;
    subplot(nr,nc,1); Nbin = plot_histo(dT_all,[-5*sig : sig/4 : 5*sig]); title('Actual data'); grid on;
    subplot(nr,nc,2); Nbin = plot_histo(sig*randn(Nmeas,1),[-5*sig : sig/4 : 5*sig]); grid on;
    title(['Gaussian distribution with sigma = ' sprintf('%.3f',sig)]);
    orient tall, wysiwyg
    
    %-------------------------------
    % create errors to add to data
    
    sig = 0.1; stag = '0p1';
    N = 10000;
   
    vals = sig*randn(N,1);
    
    figure; Nbin = plot_histo(vals,[-5*sig : sig/4 : 5*sig]);
    grid on; xlabel('Perturbation in measurement');
    title(['sigma = ' sprintf('%.3f',sig) ', std = ' sprintf('%.6f',std(vals)) ]);
    
    if iwrite==1
        ofile = [odir 'sigma_' stag '_pert.dat'];
        fid = fopen(ofile,'w');
        for ii = 1:length(vals)
            fprintf(fid,'%20.12e\n', vals(ii));   
        end
        fclose(fid);
        
        % check
        sigt = load(ofile);
        figure; Nbin = plot_histo(sigt,[-5*sig : sig/4 : 5*sig]);
        grid on; xlabel('Perturbation in measurement');
    end
end

%=============================================================

