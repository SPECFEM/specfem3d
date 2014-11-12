%
% wave2d_cg_poly.m
% Carl Tape, 21-Nov-2008
% 
% This program takes the output from wave2d.f90 and plots the polynomials
% used in the line-search conjugate gradient algorithm within the code.
%
% For GMT plotting, we execute wave2d_cg_figs.m after this program.
%
% calls
%   get_model_vec.m
%   cubic_min_4.m
%   genfit.m
%   theoryHyp.m
%   quad_shift.m
%   linefit.m
%   axes_expand.m, fontsize.m, colors.m, wysiwyg.m
% called by xxx

format long
format compact
close all
clear

% add path to additional matlab scripts
path([pwd '/matlab_scripts'],path);

colors;
npts = 100;
stfm = '%4.4i';

%---------------------------------------------------------
% USER INPUT

% base directory
dir0 = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_work/';
if ~exist(dir0), error([dir0 ' does not exist']); end

parms = [100 0 3];    % structure inversion: 1 source, Nfac=3
%parms = [200 0 3];    % source inversion (xs,ys,ts): single source
%parms = [300 0 3];    % joint inversion: 1 source, Nfac=3
%parms = [400 0 3];    % joint inversion: 5 sources, Nfac=3
%parms = [500 0 3];    % joint inversion: 25 sources, Nfac=3
%parms = [600 0 3];    % structure inversion: 25 sources, Nfac=2

parms = [6000 0 3];

%---------------------------------------------------------

irun0 = parms(1);       % irun for m0
%niter = parms(2);       % number of iterations
icubic = parms(2);      % cubic (=1) or quadratic (=0) fitting
ifit = parms(3);        % function to fit chi(m): order of polynomial

% bottom level for test parabola (=0 if no model norm term or data errors)
chi_data_stop = 0.5;

%---------------------------------------------------------

odir0 = 'OUTPUT/';
odir  = [dir0 odir0];
%odir  = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_OUTPUT_OLD_2/';  % CHT
strun0  = ['irun0 = ' num2str(irun0) ];

% determine the number of iterations
k = 0;
while k < 100
    irun = irun0 + k;
    chifile = [odir 'run_' sprintf(stfm,irun) '/chi.dat'];
    iexist = exist(chifile);
    if iexist == 2
        chi = load(chifile);
        chis_all(k+1) = chi;
    else
        break
    end
    k = k + 1;
end
istop = k - 1;
niter = ceil(istop/2);

% default is that the final model IS a test model
ifinal_test = 0;
irun_vec = [irun0+1 : 2 : irun0+2*niter];

% if the simulation converged on a test model, then tack on a repeated
% misfit value
if mod(istop,2)==1
    chis_all(istop+2) = chis_all(istop+1);
    ifinal_test = 1;
end

chis  = chis_all([1 : 2 : 2*niter+1])';
chits = chis_all([2 : 2 : 2*niter])';

% irun_vec = [irun0+1:2:irun0+2*niter];
% 
% chis = zeros(niter+1,1);
% chits = zeros(niter,1);
% 
% % load initial chi value
% stf = [odir 'run_' sprintf(stfm,irun0) '/'];
% ww = 'summed_chi_all';
% load([stf ww '.dat']); chi = eval(ww); chis(1) = chi;
% 
% % load the remaining chi values
% for ii=1:niter
%     irun = irun_vec(ii);
%     
%     stft = [odir 'run_' num2str(sprintf(stfm,irun)) '/'];
%     stf  = [odir 'run_' num2str(sprintf(stfm,irun+1)) '/'];
%     chit = load([stft 'summed_chi_all.dat']);
%     chi  = load([stf 'summed_chi_all.dat']);
%     chis(ii+1) = chi;
%     chits(ii)  = chit;
% end

%disp('         chis               chis_test');
%disp([chis chits]);
disp('chis :'); disp([chis]);
disp('chis_test :'); disp([chits]);

figure;
nc=3; nr=max( ceil(niter/nc), 2);
%nr=2; nc=2;
stniter = ['niter = ' num2str(niter) ];
stirun  = [strun0 '; ' stniter];

for ii=1:niter
    disp('----------------')
    irun = irun_vec(ii);
    stk  = num2str(ii-1);

    % load polynomial values
    % wave2d: xx1,xx2,yy1,yy2,g1,g2,Pa,Pb,Pc,Pd,xmin
    subplot(nr,nc,ii);
    stlabs = {'\nu',['S^{' stk '} ( \nu )'],'0',['\nu_{' stk '}'],['\nu_{' stk 't}']};
    %if or(ifinal_test==0, ii < niter)
    %if or(ifinal_test==1, ii < niter)
        if icubic == 1
            ww = 'cubic_poly';
            stf = [odir 'run_' num2str(sprintf(stfm,irun)) '/'];
            load([stf ww '.dat']); temp = eval(ww)';

            x1 = temp(1); x2 = temp(2);
            y1 = temp(3); y2 = temp(4);
            g1 = temp(5); g2 = temp(6);
        else
            ww = 'quad_poly';
            stf = [odir 'run_' num2str(sprintf(stfm,irun)) '/'];
            load([stf ww '.dat']); temp = eval(ww)';

            x1 = temp(1); x2 = temp(2);
            y1 = temp(3); y2 = temp(4);
            g1 = temp(5);
        end

        if icubic == 1
            [lam0,P1] = cubic_min_4(x1,x2,y1,y2,g1,g2,[1 1],stlabs);

            % check the Fortran output
            temp(7:10) ./ P1 - 1
            temp(11)/lam0 - 1
        else
            [lam0,P1] = quad_min_4(x1,x2,y1,y2,g1,[1 1],stlabs);

            % check the Fortran output
            temp(6:8) ./ P1 - 1
            temp(9)/lam0 - 1
        end
    %end
    
    % load actual chi-value computed from the next run
    % wave2d: xx1,xx2,yy1,yy2,g1,g2,Pa,Pb,Pc,Pd,xmin
    %ww = 'summed_chi_all';
    %stf = [dir0  odir 'run_' num2str(sprintf(stfm,irun+1)) '/'];
    %load([stf ww '.dat']);
    %chi = eval(ww);
    %chis(ii+1) = chi;
    
    %if ii < niter
        chi = chis(ii+1);
        msize = 6;
        axi = xlim;
        plot([lam0 lam0 axi(1)],[0 chi chi],'k');
        plot(lam0,chi,'bo','markersize',msize,'MarkerFaceColor','b');
        %axis tight; ax1 = axis; axis([ax1(1:2) 0 ax1(4)]);
        if lam0 > axi(2), axis tight; end
    %end

end
title(stirun);
disp(' chis:'); disp(chis);

if 0==1
    % figure for paper (irun0 = 20, 100)
    %save('chis','chis_quad','chis_cubic','nchi');
    load('chis');
    its = [0:nchi-1]';
    
    hess_ray = load('/home/carltape/wave2d/2d_adjoint_banana/OUTPUT_banana/run_0023/summed_chi_all.dat');
    hess_ker = load('/home/carltape/wave2d/2d_adjoint_banana/OUTPUT_banana/run_0022/summed_chi_all.dat');
    
    figure; ax1 = [-1 nchi 10^-1 10^4];
    semilogy(its,chis_quad,'r.',its,chis_cubic,'b.',...
        its,hess_ray*ones(nchi,1),'g',its,hess_ker*ones(nchi,1),'k',...
        'markersize',18,'linewidth',2);
    legend(' fit with quadratic',' fit with cubic',' hessian ray',' hessian kernel');
    axis(ax1); grid on;
    xlabel(' model number (iteration)'); ylabel(' chi(m)');
    orient tall, wysiwyg
    
    %-------------------------

    % write information to file for GMT plotting
    dir = [dir0 'gji_paper/figures/'];
    ww = 'chi_cubic_quad.dat';
    fid = fopen([dir ww],'w');
    for ii=1:nchi
        fprintf(fid,'%16.7e%16.7e%16.7e%16.7e%16.7e\n',...
            its(ii),chis_quad(ii),chis_cubic(ii),hess_ray,hess_ker);   
    end
    fclose(fid);
    
    % write information to file for GMT plotting
    ww = 'chi_cubic_quad_axes.dat';
    fid = fopen([dir ww],'w');
    fprintf(fid,'%16.7e%16.7e%16.7e%16.7e\n',ax1(1),ax1(2),ax1(3),ax1(4));   
end

if 0==1
    % 26-July-2006, comparison of convergence curves (irun0 = 4000, 4050, 4100)
    % save('chis2','chis_tt','chis_amp','chis_wav','nchi');
    %
    % 26-Oct-2006, comparison of convergence curves
    % save('chis3','chis_wav','chis_tt_xcor','chis_lnA_xcor','chis_tt_mtm','nchi');
    load('chis3');
    its = [0:nchi-1]';
    
    chis_wav = chis_wav / chis_wav(1);
    chis_tt_xcor = chis_tt_xcor / chis_tt_xcor(1);
    chis_lnA_xcor = chis_lnA_xcor / chis_lnA_xcor(1);
    chis_tt_mtm = chis_tt_mtm / chis_tt_mtm(1);
    
    %hess_ray = load('/home/carltape/wave2d/2d_adjoint_banana/OUTPUT_banana/run_0023/summed_chi_all.dat');
    %hess_ker = load('/home/carltape/wave2d/2d_adjoint_banana/OUTPUT_banana/run_0022/summed_chi_all.dat');
    
    figure; ax1 = [-1 nchi 10^-3 10^0];
    %semilogy(its,chis_tt_xcor,'r.',its,chis_lnA_xcor,'b.',its,chis_wav,'k.','markersize',18);
    %legend(' tt-xcorr (4000)',' amplitude (4050)',' waveform (4100)');
    semilogy(its,chis_wav,'k.',its,chis_tt_xcor,'r.',its,chis_lnA_xcor,'b.',its,chis_tt_mtm,'r+','markersize',18);
    legend(' waveform (4300)',' tt-xcorr (4350)',' amplitude (4400)');
    axis(ax1); grid on;
    xlabel(' model number (iteration)'); ylabel(' chi(m)');
    orient tall, wysiwyg
end

%---------------------------------------
% plot misfit as a function of iteration

nchi = length(chis);
its = [0:nchi-1]';

% x-y data to fit
x = its;
y = log10(chis);

% threshold points to remove Inf or NaN
igood = ~sum([isnan(y) isinf(y)],2);
x = x(igood);
y = y(igood);
its = its(igood);
chis = chis(igood);
nchi = length(chis);

% choose type of fit to log(chi) data: line(1), parabola(2), hyperbola(3)
%ifit = 1;
stlabs = {'line','parabola','hyperbola'};

ifac1 = 0;
ifac2 = 3;     % number of extra iterations to extrapolate
its_smooth = linspace(0-ifac1,nchi-1+ifac2,100);
sfm = '%.4e';
if ifit==1
    % fit to line for log10(chi)-its space
    %[xf,yf,mf,bf,rms] = linefit(its(2:end), log10(chis(2:end)));
    [xf,yf,mf,bf,rms] = linefit(x,y);
    yfit = mf*x + bf;
    yfit_smooth = mf * its_smooth + bf;
    
    stit1 = [' y = ( -' num2str(sprintf(sfm,mf)) ' x + ' num2str(sprintf(sfm,bf)) ')'];
    stit2 = ' ';  
    
elseif ifit==2
    % fit to parabola for log10(chi)-its space
    [xf,yf,P,rms,stit] = parabolafit(x,y);
    yfit = P(1)*x.^2 + P(2)*x + P(3);
    yfit_smooth = P(1)*its_smooth.^2 + P(2)*its_smooth + P(3);
    
elseif ifit==3
    % initial guesses based on a line
    [xf,yf,mf,bf,rms] = linefit(x,y);
    m0 = [-mf 0 1 -bf]';
    y0 = theoryHyp(m0, x);
    figure; hold on; plot(x,y,'.',x,y0,'ro');
    
    % perturbations
    jogvec = [1e-8 * ones(4,1)]';

    disp('  ');
    disp('Model: A,B,C,D in equation Ax + Bxy + Cy + D = 0');
    disp('                   ---> y = (-Ax - D) / (Bx + C)');
    mz = m0; itmx = 5;
    for ii=1:itmx
        [mz e1] = genfit('theoryHyp', mz, jogvec, y, [x]);
        mz
        yest = theoryHyp(mz,x);
        res = y - yest; rms = sqrt( (res' * res) / length(res) );
        stRMS = [' RMS = ' num2str(sprintf('%.4e', rms)) ';'];
    end
    disp('Best-fit hyperbola:');
    mz
    a = mz(1); b = mz(2); c = mz(3); d = mz(4); 
    stit1 = [' y = ( -' num2str(sprintf(sfm,a)) ' x - ' num2str(sprintf(sfm,d)) ') / (' ...
                       num2str(sprintf(sfm,b)) ' x + ' num2str(sprintf(sfm,c)) ')'];
    stit2 = [' y = ' num2str(sprintf(sfm,a)) ' x + ' num2str(sprintf(sfm,b)) ' xy + ' ...
            num2str(sprintf(sfm,c)) ' y + ' num2str(sprintf(sfm,d)) ' = 0'];               
    disp(stRMS);
    
    yfit = theoryHyp(mz,x);
    yfit_smooth = theoryHyp(mz,its_smooth);
    plot(its_smooth, yfit_smooth, 'r--');
    legend(' data',' initial guess line',' hyperbola fit');
end
chifit = 10.^yfit;
chifit_smooth = 10.^yfit_smooth;

%------------------------------

figure; nr=2; nc=2;
msize = 24; lsize = 2;
xlims = [-1 its_smooth(end)];
xticks = [round(its_smooth(1)):round(its_smooth(end))];
stit = {[' fitting a ' stlabs{ifit} ' to log10(chi)-vs-iteration data'],stit1,stit2};

stx = ' k, model number (iteration)';

subplot(nr,nc,1); hold on;
plot(its_smooth,chifit_smooth,'r--','linewidth',lsize);
plot(its,chis,'.','markersize',msize);
xlim(xlims); grid on; set(gca,'xtick',xticks);
xlabel(stx);
ylabel(' S (m)','fontsize',18);
title(stirun);

subplot(nr,nc,2);
if 1==1
    plot(its_smooth,log10(chifit_smooth),'r--',its,log10(chis),'.','markersize',msize,'linewidth',lsize);
    ylabel(' log10 [ S (m) ]','fontsize',18);
    %set(gca,'ytick',[-10:10]);
else
    semilogy(its_smooth,chifit_smooth,'r--',its,chis,'.','markersize',msize,'linewidth',lsize);
    ylabel(' S (m)','fontsize',18);
    set(gca,'ytick',10.^[-10:10]);
end
grid on; xlim(xlims); set(gca,'xtick',xticks);
xlabel(stx); title(stit);

% variance reduction
chi_before = chis(1:nchi-1);
chi_after  = chis(2:nchi);
%var_red1   = 100 * ( 1 - ( (chi_after-chi_data_stop).^2 ./ (chi_before-chi_data_stop).^2 ) );
%var_red2   = 100 * ( 1 - ( (chi_after-chi_data_stop).^2 ./ (chi_before-chi_data_stop).^2 ) );
var_red1   = log( chi_before ./ chi_after );
var_red2   = log( chi_before ./ chi_after );

disp('  '); disp('VARIANCE REDUCTION');
disp([its(2:end) chi_before chi_after var_red1]);

% var red corresponding to best-fitting chi-vs-its curve
its2         = xticks;
chifit2      = interp1(its_smooth,chifit_smooth,its2);
%var_red1_fit = 100 * ( 1 - ( chifit2(2:end).^2 ./ chifit2(1:end-1).^2 ) );
var_red1_fit = log( chifit2(1:end-1) ./ chifit2(2:end) )

x_var     = its(2:end);
x_var_fit = its2(2:end);

subplot(nr,nc,3); hold on;
plot(x_var_fit,var_red1_fit,'r.--','markersize',msize,'linewidth',lsize);
plot(x_var,var_red1,'b.-','markersize',msize,'linewidth',lsize);
axis([xlims 0 1]); grid on; set(gca,'xtick',xticks);
xlabel(stx);
ylabel(' Logarithmic Variance Reduction');
title(' Variance reduction between successive models');

% subplot(nr,nc,4); plot(its(2:end),var_red2,'.-');
% xlim(xlims); grid on; set(gca,'xtick',its);
% xlabel(stx);
% ylabel(' variance reduction');

lchi_smooth = log10(chifit_smooth);
converge_fit = abs( gradient( lchi_smooth, its_smooth) );
converge_its = 0.5+its(1:end-1);
converge_pts = abs( diff(log10(chis)) );

subplot(nr,nc,4);
%semilogy(its_smooth,-gradient(chifit_smooth, its_smooth),'r--','linewidth',lsize);
plot(its_smooth,converge_fit,'r--',converge_its,converge_pts,'.','markersize',msize,'linewidth',lsize);
axis([xlims 0 0.7]); grid on; set(gca,'xtick',xticks);
xlabel(stx);
ylabel(' Order of convergence rate, | \Delta log10(S) /  \Delta k |');

fontsize(9), orient tall, wysiwyg

%-----------------------------
% This generates the initial polynomial, in order to have the values
% of the fitting curves to write to file in gji_figs.m for GMT plotting.

% specify number of iterations to plot; load initial chi
niter = 1;
stirun = [' irun0 = ' num2str(irun0) '; niter = ' num2str(niter) ];
stf = [odir 'run_' num2str(sprintf(stfm,irun0)) '/'];
ww = 'chi';
load([stf ww '.dat']); chi = eval(ww); chis(1) = chi;

%figure; %nr=3; nc=3;
for ii=1:niter
    irun = irun_vec(ii);
    stk  = num2str(ii-1);
    stlabs = {'\nu',['S^{' stk '} ( \nu )'],'0',['\nu_{' stk '}'],['\nu_{' stk 't}']};
    
    if icubic == 1
        ww = 'cubic_poly';
        stf = [odir 'run_' num2str(sprintf(stfm,irun)) '/'];
        load([stf ww '.dat']); temp = eval(ww)';
        x1 = temp(1); x2 = temp(2);
        y1 = temp(3); y2 = temp(4);
        g1 = temp(5); g2 = temp(6);
    else
        ww = 'quad_poly';
        stf = [odir 'run_' num2str(sprintf(stfm,irun)) '/'];
        load([stf ww '.dat']); temp = eval(ww)';
        x1 = temp(1); x2 = temp(2);
        y1 = temp(3); y2 = temp(4);
        g1 = temp(5);
    end
    
    if 0==1     % call functions for plotting
        if icubic == 1
            [xmin,P1] = cubic_min_4(x1,x2,y1,y2,g1,g2,[1 0],stlabs);
        else
            [xmin,P1] = quad_min_4(x1,x2,y1,y2,g1,[1 0],stlabs);
        end
        
    else        % call functions as scripts
        
        % quadratic polynomial for computing test-model nu
        aquad = g1^2/(4*y1);
        Pquad_test = [aquad g1 y1]';
        n = 100;
        
        % quadratic interpolation (only used to get chi if icubic = 1)
        a = ((y2 - y1) - g1*(x2 - x1)) / (x2^2 - x1^2);
        b = g1;
        c = y1 - a*x1^2 - b*x1;
        Pquad = [a b c]';
        [Pquad2,qvert,stit] = quad_shift(Pquad,1);
        
        if icubic == 1  % copied from cubic_min_4.m on 2-14-06
            
            a = ( -2*(y2-y1) + (g1+g2)*(x2-x1) ) / (x2 - x1)^3;
            b = ( 3*(y2-y1) - (2*g1 + g2)*(x2 - x1) ) / (x2 - x1)^2;
            c = g1;
            d = y1;
            P2 = [a b c d]';
            Pcubic = cubic_shift(P2,x1,0);
            xmin = cubic_min(Pcubic,x1);
            
            figure; hold on;
            specs = [2 14 18 12]; fac = 0.05;
            axpoly = axes_expand([x1 x2 0 max([y1 y2])],1.2,1);
            axpoly(3) = 0;
            dy = axpoly(4) - axpoly(3);
            dx = axpoly(2) - axpoly(1);
            ylab = axpoly(3) - fac*dy;
            ymin = polyval(Pcubic,xmin);
            
            xpts = linspace(axpoly(1),axpoly(2),n);    
            g1_line = polyval([g1 y1-g1*x1],xpts);
            g2_line = polyval([g2 y2-g2*x2],xpts);
            g1_quad_test = polyval(Pquad_test,xpts);
            g1_quad_fit = polyval(Pquad,xpts);
            g1_cube_fit = polyval(Pcubic,xpts);
            
            plot(xpts,g1_cube_fit,'b','linewidth',specs(1));
            plot(xpts,g1_quad_fit,'r--','linewidth',specs(1));
            plot(xpts,g1_quad_test,'b--');
            plot(xpts,g1_line,'r--',xpts,g2_line,'r--');
            
            plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
            plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
            plot([xmin xmin axpoly(1)],[0 ymin ymin],'k--');
            
            plot([x1 x2],[y1 y2],'ko','markersize',specs(2),'MarkerFaceColor','b');
            plot(xmin,ymin,'bo',x2,0,'bo',xmin,0,'bo','markersize',specs(2),'MarkerFaceColor','w');
            plot(0.5*x2,0,'go','markersize',specs(2),'MarkerFaceColor','w');  % linear projection
            
            axis(axpoly);
            xlabel(stlabs{1},'fontsize',specs(3));
            ylabel(stlabs{2},'fontsize',specs(3));
            grid on;
            text(xmin,ylab,stlabs{4},'fontsize',specs(4));
            text(x2,ylab,stlabs{5},'fontsize',specs(4));
            orient tall, wysiwyg
            
        else  % copied from quad_min_4.m on 2-14-06
            
            xmin = qvert(1);

            figure; hold on;
            specs = [2 14 18 12]; fac = 0.05;
            axpoly = axes_expand([x1 x2 0 max([y1 y2])],1.2,1);
            axpoly(3) = 0;
            dy = axpoly(4) - axpoly(3);
            dx = axpoly(2) - axpoly(1);
            ylab = axpoly(3) - fac*dy;
            ymin = polyval(Pquad,xmin);
            
            xpts = linspace(axpoly(1),axpoly(2),n);
            g1_line = polyval([g1 y1-g1*x1],xpts);
            g2_line = zeros(n,1);
            g1_quad_test = polyval(Pquad_test,xpts);
            g1_quad_fit = polyval(Pquad,xpts);
            g1_cube_fit = zeros(n,1);
            
            plot(xpts,g1_quad_fit,'b','linewidth',specs(1));
            plot(xpts,g1_quad_test,'b--','linewidth',specs(1));
            plot(xpts,g1_line,'r--','linewidth',specs(1));
            
            plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
            plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
            plot([xmin xmin axpoly(1)],[0 ymin ymin],'k--');
            plot([x1 x2],[y1 y2],'ko','markersize',specs(2),'MarkerFaceColor','b');
            plot(xmin,ymin,'bo',x2,0,'bo',xmin,0,'bo','markersize',specs(2),'MarkerFaceColor','w');
            plot(0.5*x2,0,'go','markersize',specs(2),'MarkerFaceColor','w');  % linear projection
            axis(axpoly);
            xlabel(stlabs{1},'fontsize',specs(3));
            ylabel(stlabs{2},'fontsize',specs(3));
            title({stit{1},stit{2}},'fontsize',specs(3))
            grid on;
            text(xmin,ylab,stlabs{4},'fontsize',specs(4));
            text(x2,ylab,stlabs{5},'fontsize',specs(4));
            orient tall, wysiwyg
        end
    end
        
    % load actual chi-value computed from the next run
    % wave2d: xx1,xx2,yy1,yy2,g1,g2,Pa,Pb,Pc,Pd,xmin
    ww = 'chi';
    stf = [odir 'run_' num2str(sprintf(stfm,irun+1)) '/'];
    load([stf ww '.dat']);
    chi = eval(ww);
    
    msize = 14;
    axi = xlim; plot([xmin xmin axi(1)],[0 chi chi],'k');
    plot(xmin,chi,'bo','markersize',msize,'MarkerFaceColor','b');
end
title(stirun);
orient tall, wysiwyg

%=============================================================
