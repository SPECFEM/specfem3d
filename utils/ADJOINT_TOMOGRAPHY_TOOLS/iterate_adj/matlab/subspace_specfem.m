%
% subspace_specfem.m
% CARL TAPE, 01-Feb-2009
% printed xxx
%
% This program assumes you have FIRST run compute_misfit.m, and also have
% generated the corresponding Hessian matrix to load here.
%
% See also wave2d_subspace.m
% 
% calls hfile2hess.m
% called by xxx
%

clear
close all
format short
format compact

% add path to additional matlab scripts
path(path,[pwd '/matlab_scripts']);

ax1 = [-121 -114 31 37];
stfm = '%4.4i';

xlab1 = 'p, singular value index';
xlab2 = 'p, singular value truncation index';
ylab1 = 'log10( singular value )';
ylab2 = 'log10( misfit : dot[ d - G*dm(p), d - G*dm(p) ] )';

imod = input(' Enter the current model number (0, 1, ...): ');
iwrite = input(' Enter 1 to write files (0 otherwise) : ');

stmod = sprintf('m%2.2i',imod);

% directory containing the Hessian file that was output from the cluster
sthess = ['/home/carltape/results/HESSIANS/'];

% load the data and event weights (compute_misfit.m)
idatacov = input(' Enter idatacov (1--by event, 2--by window, 3--none) : ');
stdcovs = {'event','window','none'};
stdcov = stdcovs{idatacov};
ftag = [stmod '_' stdcov];
odir0 = ['OUTPUT_SUBSPACE/' stmod '/' stdcov '/'];
load([odir0 ftag '_data_norms']);
nsrc = length(eids);

%----------------------------------------------
% check the slick operation used in subspace_update.f90

if 0==1
    npar = 10;
    nsrc = 4;
    G = rand(npar,nsrc);
    d = rand(nsrc,1);
    
    x = zeros(npar,1);
    for isrc=1:nsrc      % loop over kernels (hundreds)
       for ipar=1:npar   % loop over gridpoints (millions)
          x(ipar) = x(ipar) + G(ipar,isrc) * d(isrc); 
       end
    end
    G*d, x
end

%----------------------------------------------

% construct the matrix of data-covariance normalization
Cd_fac = zeros(nsrc,nsrc);
for i = 1:nsrc
    for j = 1:nsrc
        Cd_fac(i,j) = dcov_fac_e(i) * dcov_fac_e(j);
        i_ind(i,j) = i;
        j_ind(i,j) = j;
    end
end

% % OBSOLETE: construct the matrix of weight factors and data-covariance normalization
% if 0==1
%     Ws = zeros(nsrc,nsrc);
%     Cd_fac = zeros(nsrc,nsrc);
%     for i = 1:nsrc
%         for j = 1:nsrc
%             Ws(i,j) = ws(i) * ws(j);
%             Cd_fac(i,j) = dcov_fac_e(i) * dcov_fac_e(j);
%             i_ind(i,j) = i;
%             j_ind(i,j) = j;
%         end
%     end
%     figure; pcolor(i_ind,j_ind,Ws); shading flat;
%     xlabel('Row index'); ylabel('Column index');
%     title('Symmetric matrix of data weight terms');
%     %caxis([0 200]);
%     colorbar; axis equal; axis tight; 
% end

% load the event IDs corresponding to the kernels
%eids = load([sthess 'kernels_done_mod']);
%neid = length(eids);

kmin = 6;
kmax = 6;
nker = kmax - kmin + 1;
stkers = {'mu_kernel','kappa_kernel','both','mu_kernel_smooth','kappa_kernel_smooth','both_smooth'};
stkerlabs = {'unsmoothed beta kernel','unsmoothed beta kernel','beta + bulk','smoothed beta kernel','smoothed bulk kernel','smoothed beta + bulk'};

%stwts = {'with data weights'};
%stwts = {'with data weights','with NO data weights'};
nker = length(stkers);
%stcs = {'b.','r.'};
%iw = 1;

%for iw = 1:length(stwts)

for iker = kmin:kmax
    stker = stkers{iker};
    stkerlab = stkerlabs{iker};
    stlab = [ftag '_' stker];

    dir1 = [odir0 stker '/'];
    
    if ~any(iker == [3 6]) 
        % get the Hessian from the subspace_hessian.f90 output
        hfile = sprintf([sthess 'hessian_index_' stmod '_' stker '_%3.3i'],nsrc);
        H0 = hfile2hess(hfile,nsrc);
        H = H0;

        % include the weights determined from the data
        %if iw == 1
        %    H = H .* Ws;
        %end

        % When we computed the kernels, we did not include the
        % normalization factor based on the number of windows picked
        % that is included within the data covariance (and used in computing Ws).
        % Thus we compute the MATRIX Cd_fac above to account for this.
        H = H ./ Cd_fac;

        % put the matrix onto order-1 values for computation purposes
        %Hnorm = norm(H); H = H / norm(H);
        
    else
        hfile_beta = sprintf([sthess 'hessian_index_' stmod '_' stkers{iker-2} '_%3.3i'],nsrc);
        H0_beta = hfile2hess(hfile_beta,nsrc);
        hfile_bulk = sprintf([sthess 'hessian_index_' stmod '_' stkers{iker-1} '_%3.3i'],nsrc);
        H0_bulk = hfile2hess(hfile_bulk,nsrc);
        
        %H_beta = H0_beta .* Ws ./ Cd_fac;
        %H_bulk = H0_bulk .* Ws ./ Cd_fac;
        H_beta = H0_beta ./ Cd_fac;
        H_bulk = H0_bulk ./ Cd_fac;
        H = H_beta + H_bulk;
        
        dtemp = [[1:nsrc]' eids_num diag(H_beta) diag(H_bulk) diag(H) diag(H_bulk)./diag(H_beta) Ns];
        dtemp_sort1 = sortrows(dtemp,-6);   % sort by relative importance of bulk-to-shear
        dtemp_sort2 = sortrows(dtemp,-5);   % sort by diagonal
        
        % sort all entries of the Hessian -- note that some are negative
        % we take the upper traingular elements, excluding the diagonal
        Htri = triu(H,1); itri = triu(i_ind,1); jtri = triu(j_ind,1);
        Hcols = [Htri(:) itri(:)  jtri(:) ];
        [Hsort,iHsort_good] = sortrows(Hcols, -1);
        [Hsort,iHsort_bad] = sortrows(Hcols, 1);
        
        if iwrite == 1
            filename = [dir1 'hessian_bulk-shear_' stlab];
            fid = fopen(filename,'w');
            fprintf(fid,'Hessian diagonal contributions from beta and bulk:\n');
            fprintf(fid,'%6s%10s%10s%10s%10s%10s%6s\n','index','eid','beta','bulk','total','bulk/beta','Nwin');
            for ix = 1:nsrc
                fprintf(fid,'%6i%10i%10.2e%10.2e%10.2e%10.4f%6i\n',dtemp_sort1(ix,:));
            end
            fclose(fid);
            
            filename = [dir1 'hessian_overall_' stlab];
            fid = fopen(filename,'w');
            fprintf(fid,'Hessian diagonal contributions from beta and bulk:\n');
            fprintf(fid,'%6s%10s%10s%10s%10s%10s%6s\n','index','eid','beta','bulk','total','bulk/beta','Nwin');
            for ix = 1:nsrc
                fprintf(fid,'%6i%10i%10.2e%10.2e%10.2e%10.4f%6i\n',dtemp_sort2(ix,:));
            end
            fclose(fid);
            
            % make a list of the largest N positive and negative Hessian elements
            nprint = 100;      % must be less than nsrc*nsrc
            filename = [dir1 'hessian_elements_good_' stlab];
            fid = fopen(filename,'w');
            fprintf(fid,'Largest positive non-diagonal elements of Hessian:\n');
            fprintf(fid,'%10s%6s%6s%12s%12s\n','Hij','i','j','eid-i','eid-j');
            for ix = 1:nprint
                iH = iHsort_good(ix);
                fprintf(fid,'%10.2e%6i%6i%12s%12s\n',Hcols(iH,:),eids{Hcols(iH,2)},eids{Hcols(iH,3)});
            end
            fclose(fid);
            
            filename = [dir1 'hessian_elements_bad_' stlab];
            fid = fopen(filename,'w');
            fprintf(fid,'Largest negative non-diagonal elements of Hessian:\n');
            fprintf(fid,'%10s%6s%6s%12s%12s\n','Hij','i','j','eid-i','eid-j');
            for ix = 1:nprint
                iH = iHsort_bad(ix);
                fprintf(fid,'%10.2e%6i%6i%12s%12s\n',Hcols(iH,:),eids{Hcols(iH,2)},eids{Hcols(iH,3)});
            end
            fclose(fid);
        end
        
        disp(' gradient balance for the Hessian (subspace) inversion (mean of last column) :');
        disp(mean( diag(H_bulk)./diag(H_beta) ))
    end
    
    %H = H + eye(nsrc);
    % mutemp = inv(H + eye(nsrc)) * dnorm;

    dH = diag(H);
    trH = sum(dH);
    trHH = sum(diag(H'*H));
    
    %------------------------
    
    pinds = [1:nsrc]';
    
    disp(' properties of Hessian (min, median, mean(abs), max, std):');
    stH1 = sprintf('min %.1e, median %.1e, mean(abs) %.1e, max %.1e, std %.1e, sum %.1e',...
        min(H(:)), median(H(:)), mean(abs(H(:))), max(H(:)), std(H(:)), sum(H(:)) );
    stH2 = sprintf('DIAGONAL : min %.1e, median %.1e, mean(abs) %.1e, max %.1e, std %.1e, sum %.1e',...
        min(dH), median(dH), mean(abs(dH)), max(dH), std(dH), sum(dH) );
    disp(stH1);
    disp(stH2);
    
    % plot H
    Hfac = 0.5;    % adjust to saturate the color scale as you wish
    figure; pcolor(i_ind,j_ind,H); shading flat;
    xlabel('Row index'); ylabel('Column index');
    title({['Hessian (symmetric matrix) for ' stkerlab],stH1,stH2});
    caxis([-1 1]*max(H(:))*Hfac), colormap('jet'), colorbar;
    %map = colormap('seis'); colormap(flipud(map));
    axis equal; axis tight; 
    fontsize(10), orient tall, wysiwyg

    % sort the diagonal values of H
    Hdiag = diag(H);
    %Hdiag = diag(H0);
    [jk,isort] = sortrows([pinds eids_num Hdiag],-3);
    npick = min([20 nsrc]);
    npick = nsrc;
    disp('  '); disp([' first ' num2str(npick) ' sorted diagonal entries of H: ' stkerlab]);
    for ii = 1:npick
        jj = isort(ii);
       disp(sprintf('%3i %10i (%3i), Ne = %5i, Hkk  =  %10.4e',...
           ii,eids_num(jj),pinds(jj),Ns(jj),Hdiag(jj)));
    end
    
    [jtemp,isort] = sort(Hdiag,'descend');
    
    figure(20+iker); hold on;
    plot(pinds,log10(sort(Hdiag,'descend')),'ro');
    %plot(pinds,log10(sort(diag(H'*H),'descend')),'b.');
    plot([1 nsrc],log10(trH)*[1 1],'r--');
    %plot([1 nsrc],log10(trHH)*[1 1],'b--');
    %legend('diagonal element of H','diagonal element of HtH','trace (H)','trace (HtH)');
    legend('diagonal element of H','trace (H)');
    xlabel('sorted event index'); ylabel(' log10 H-kk');
    title(sprintf(' trace of H = %6.2e --- trace of HtH = %6.2e',trH,trHH));
    xlim([1 nsrc]); grid on;
    fontsize(11), orient tall, wysiwyg
    
    %----------------------------------------------
    
    itsvd = 0;

    if itsvd == 0

        % regularization choices
        numlam = 100;
        if 0==1
            minlampwr = -2.25; maxlampwr = 1;
            lampwr = linspace(minlampwr,maxlampwr,numlam);
            lamvec = 10.^lampwr * sum(dH);                  % based on trace of H
            %lamvec = 10.^lampwr * sqrt(sum(diag(H'*H)));    % based on trace of H'*H
        else
            %minlampwr = -2; maxlampwr = 5;                 % m03
            minlampwr = -3; maxlampwr = 3;                 % m04 - m11
            minlampwr = -5; maxlampwr = 0;                 % m12
            lampwr = linspace(minlampwr,maxlampwr,numlam);
            lamvec = 10.^lampwr;
        end

        [f_h, rss, mss, Gvec, Fvec, dof, kap, iL, iGCV, iOCV] = ridge_carl(dnorm,H,lamvec);
        title({['model ' stmod],stkerlab,sprintf('iL = %i,  iOCV = %i,  iGCV = %i',iL,iOCV,iGCV)});
        
        mu_all = zeros(nsrc,numlam);
        mu_norm = zeros(numlam,1);
        mu_res = zeros(numlam,1);
        for ip = 1:numlam
            lam = lamvec(ip);
            %if ip==iGCV, disp('iGCV -- GCV lambda pick'); end
            %if ip==iOCV, disp('iOCV -- OCV lambda pick'); end
            %if ip==iL, disp('iL -- L-curve lambda pick'); end
            disp(sprintf('%i out of %i, lam = %.2e',ip,numlam,lam));
                        
            mu_all(:,ip) = inv(H + lam*eye(nsrc))*dnorm;
            %mu_all(:,ip) = inv(H'*H + lam^2*eye(nsrc))*H'*dnorm;
            
            % f  = inv(G'*diag(Wvec)*G + lam0^2*eye(ngrid*ndim))*G'*diag(Wvec)*d;
            mu_norm(ip) = norm(mu_all(:,ip));
            mu_res(ip) = norm( dnorm - H*mu_all(:,ip) );
        end
        
        ipick = iOCV;     % KEY: select a damping parameter (iL, iOCV, iGCV)
        ltemp = lamvec(ipick)
        mtemp = mu_all(:,ipick);
        mtemp_norm = mu_norm(ipick)
        
        figure; nr=2; nc=2;
        
        subplot(nr,nc,1); hold on; grid on;
        plot( log10(lamvec), log10(mu_norm),'b.','markersize',10);
        plot( log10(lamvec(ipick)), log10(mu_norm(ipick)),'ro','markersize',10);
        xlabel('log10( lambda )'); ylabel('log10( norm of mu )');
        
        subplot(nr,nc,2); hold on; grid on;
        plot( log10( lamvec ), mu_res, '.');
        plot( log10( lamvec(ipick) ), mu_res(ipick),'ro','markersize',10);
        xlabel('log10( lambda )'); ylabel('dnorm residual');
        
        subplot(nr,nc,3); plot( dnorm, H*mtemp, '.');
        xlabel('dnorm obs'); ylabel('dnorm pred'); grid on; axis equal;
        
        subplot(nr,nc,4); stdres = std(dnorm - H*mtemp);
        plot_histo(dnorm - H*mtemp,[-4:0.5:4]*stdres); 
        %plot( 100*(dnorm - H*mtemp)./dnorm ,'.');
        %ylabel('dnorm residual = 100(pred - obs)/obs');
        title(' residuals :  d - H*mu');
        xlabel('event index'); grid on;

        orient tall, wysiwyg, fontsize(10)
        
%                 iCV = input(' Enter 1 for OCV or 0 for GCV : ');
%                 if iCV == 1
%                     lam = lamvec(iOCV);
%                 elseif iCV == 0
%                     lam = lamvec(iGCV);
%                 else
%                     error(' iCV must be 1 or 0');
%                 end
%                 mu = inv(H'*H + lam^2*eye(nsrc,nsrc))*H'*dnorm;   % KEY vector
%                 stp = sprintf('lambda = %.3f',lam);

         % write mu vectors to file
        if iwrite == 1
            otag = ['mu_' stlab];
            odir = [dir1 'mu_all/'];
            for ip = 1:numlam
                lam = lamvec(ip);
                mutemp = inv(H'*H + lam^2*eye(nsrc,nsrc))*H'*dnorm;
                filename = sprintf([otag '_p%3.3i'],ip);
                fid = fopen([odir filename],'w');
                for ii = 1:nsrc
                    fprintf(fid,'%24.12e\n',mutemp(ii));
                end
                fclose(fid);
            end

            fid = fopen([odir '/lambda_' stdcov],'w');
            for ip = 1:numlam
                lam = lamvec(ip);
                fprintf(fid,'%24.12e%24.12e%6i\n',lamvec(ip),1/lamvec(ip),ip);
            end
            fclose(fid);
        end

    else

        % analyze the singular values of H
        [U,S,V] = svd(H);
        s = diag(S);

        pmin = min([13 nsrc]);
        ifit = [(pmin+1):60];
        %ifit = pinds;
        [xf,yf,mf,bf,rms] = linefit(pinds(ifit), log10(s(ifit)));

        %subplot(nr,nc,iker);
        figure; hold on;
        plot(xf,yf,'g','linewidth',6);
        plot(pinds,log10(s),'b.','markersize',10);
        plot(pmin*[1 1],[min(log10(s)) max(log10(s))],'k','linewidth',2)
        grid on; xlabel(xlab1); ylabel(ylab1);
        title(['singular value spectra for Hessian : ' stker]);
        fontsize(10), orient tall, wysiwyg

        %-------------------------

        %tpinds = [1:10:100];
        tpinds = pinds;
        [mu_all,rss,f_r_ss] = tsvd(dnorm,H,tpinds);

        % norms of mu vectors
        mu_norm = zeros(nsrc,1);
        for ip = 1:nsrc
            mu_norm(ip) = norm(mu_all(:,ip));
        end

        figure;
        subplot(2,2,1); hold on; grid on;
        plot(pinds,log10(s),'b.','markersize',10);
        xlabel(xlab1); ylabel(ylab1);
        subplot(2,2,2); hold on; grid on;
        plot( tpinds, log10(rss),'b.','markersize',10)
        xlabel(xlab2); ylabel(ylab2);
        subplot(2,2,3); hold on; grid on;
        plot( tpinds, log10(mu_norm),'b.','markersize',10)
        xlabel(xlab2); ylabel('log10( norm of mu )');
        fontsize(10), orient tall, wysiwyg

        % write mu vectors to file
        if iwrite == 1
            odir = [dir1 'mu_all/'];
            for ip = 1:length(tpinds);
                pmax = tpinds(ip);
                filename = sprintf([odir '_p%3.3i'],pmax);
                fid = fopen([odir '/' filename],'w');
                for ii = 1:nsrc
                    fprintf(fid,'%24.12e\n',mu_all(ii,ip));
                end
                fclose(fid);
            end
        end
    end
    
    
end

%end

% write misfit values to file -- these do not depend on the choice of model
% variable used in constructing the Hessian
if iwrite == 1
    % write data-norm vector to file
    fid = fopen([odir0 ftag '_dmisfit'],'w');
    fprintf(fid,'%24.12e\n',dmisfit);
    fclose(fid);
    
    % write data-norm vector to file
    fid = fopen([odir0 ftag '_nwin_tot'],'w');
    fprintf(fid,'%.0f\n',N);
    fclose(fid);
    
    % write data-norm vector to file
    fid = fopen([odir0 ftag '_data_norm'],'w');
    for isrc = 1:nsrc
        fprintf(fid,'%24.12e\n',dnorm(isrc));
    end
    fclose(fid);

    % write the factors for the data covariance (nevent x 1)
    % WRITE THIS IN DOUBLE PRECISION, JUST FOR SIMPLICITY
    fid = fopen([odir0 ftag '_dcov_fac'],'w');
    for isrc = 1:nsrc
        fprintf(fid,'%24.12e\n',dcov_fac_e(isrc));
    end
    fclose(fid);
end

%======================================================================
