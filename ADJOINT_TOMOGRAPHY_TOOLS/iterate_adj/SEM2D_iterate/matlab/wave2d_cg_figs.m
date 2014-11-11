%
% wave2d_cg_figs.m
% Carl Tape, 21-Nov-2008
%
% This program is called following wave2d_cg_poly.m.
%
% calls xxx
% called by xxx

%irun0 = 140;

stfm = '%4.4i';

% axes limits
n1 = 0;
n2 = 8;
n2 = max([17 xlims(2)]);
%ax_chi = [n1-1 n2 10^0 axpoly(4)];
%ax_chi = [n1-1 n2 10.^[0 4]];
ax_chi = [n1-1 n2 10.^[-1 3] ];
ax_var = [n1-1 n2 0 100];

%n2 = 4; ax_chi = [n1-1 n2 10^-2 10^3];

% dir0 from wave2d_cg_poly.m
dir1 = [dir0 '/PLOTTING/DATA_FILES/'];
if ~exist(dir1), error([dir1 ' does not exist']); end

%==================================================================

ww = ['poly_run' num2str(sprintf(stfm,irun0))];
if 1==1
    save(ww,'ax_chi','ax_var','its','chis','its_smooth','chifit_smooth',...
        'x_var_fit','var_red1_fit','x_var','var_red1',...
        'axpoly','xpts','g1_line','g1_quad_test','g1_cube_fit','g2_line','g1_quad_fit',...
        'x1','y1','x2','y2','xmin','ymin','chi');
    %break
else
    load(ww);
end

specs = [1 6 12 10];
xticks = [n1:n2];

if 1==1    
    % write curves to file
    ww = ['poly_curve_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    for ii=1:length(xpts)
        fprintf(fid,'%16.7e%16.7e%16.7e%16.7e%16.7e%16.7e\n',...
            xpts(ii),g1_line(ii),g2_line(ii),g1_quad_test(ii),g1_cube_fit(ii),g1_quad_fit(ii) );   
    end
    fclose(fid);
    
    % write points to file
    ww = ['poly_points_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    xptsQ = [x1 x2 x2 xmin xmin xmin   0];
    yptsQ = [y1  0 y2    0 ymin  chi chi];
    for ii=1:length(xptsQ)
        fprintf(fid,'%16.7e %16.7e \n',xptsQ(ii),yptsQ(ii) );   
    end
    fclose(fid);
    
    % write chi fit to file
    ww = ['chi_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    for ii=1:length(chis)
        fprintf(fid,'%16.7e %16.7e \n',its(ii),chis(ii) );   
    end
    fclose(fid);
    
    % write chi fit to file
    ww = ['chi_curve_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    for ii=1:length(chifit_smooth)
        fprintf(fid,'%16.7e %16.7e \n',its_smooth(ii),chifit_smooth(ii) );   
    end
    fclose(fid);
    
    % write var fit to file
    ww = ['var_fit_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    for ii=1:length(var_red1_fit)
        fprintf(fid,'%16.7e %16.7e \n',x_var_fit(ii),var_red1_fit(ii) );   
    end
    fclose(fid);
    
    % write var to file
    ww = ['var_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    for ii=1:length(var_red1)
        fprintf(fid,'%16.7e %16.7e \n',x_var(ii),var_red1(ii) );   
    end
    fclose(fid);
    
    % write axes to file
    ww = ['axes_' num2str(sprintf(stfm, irun0)) '.dat'];
    fid = fopen([dir1 ww],'w');
    fprintf(fid,'%16.7e %16.7e %16.7e %16.7e \n',axpoly(1),axpoly(2),axpoly(3),axpoly(4));
    fprintf(fid,'%16.7e %16.7e %16.7e %16.7e \n',ax_chi(1),ax_chi(2),ax_chi(3),ax_chi(4));
    fprintf(fid,'%16.7e %16.7e %16.7e %16.7e \n',ax_var(1),ax_var(2),ax_var(3),ax_var(4));
    fclose(fid);
end
    
figure; nr=5; nc=4;
subplot(nr,nc,1);  hold on; 

subplot(nr,nc,2);  hold on; 

subplot(nr,nc,3);  hold on;

subplot(nr,nc,4);  hold on; axis(axpoly);
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot(x1,y1,'ko','markersize',specs(2),'MarkerFaceColor','b');

subplot(nr,nc,5);  hold on; 

subplot(nr,nc,6);  hold on; axis(axpoly);
plot(xpts,g1_line,'r--');
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot(x1,y1,'ko','markersize',specs(2),'MarkerFaceColor','b');

subplot(nr,nc,7);  hold on; axis(axpoly);
plot(xpts,g1_line,'r--');
plot(xpts,g1_quad_test,'b--');
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot(x1,y1,'ko','markersize',specs(2),'MarkerFaceColor','b');
plot(x2,0,'bo','markersize',specs(2),'MarkerFaceColor','w');

subplot(nr,nc,8);  hold on; 

subplot(nr,nc,9);  hold on; 

subplot(nr,nc,10); hold on; axis(axpoly);
plot(xpts,g1_line,'r--');
plot(xpts,g1_quad_test,'b--');
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
plot(x1,y1,'ko','markersize',specs(2),'MarkerFaceColor','b');
plot(x2,0,'bo','markersize',specs(2),'MarkerFaceColor','w');
plot(x2,y2,'ko','markersize',specs(2),'MarkerFaceColor','b');

subplot(nr,nc,11); hold on; 

subplot(nr,nc,12); hold on; axis(axpoly);
plot(xpts,g1_line,'r--',xpts,g2_line,'r--');
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
plot([x1 x2],[y1 y2],'ko','markersize',specs(2),'MarkerFaceColor','b');
plot(x2,0,'bo','markersize',specs(2),'MarkerFaceColor','w');

subplot(nr,nc,13); hold on; axis(axpoly);
plot(xpts,g1_cube_fit,'b','linewidth',specs(1));
plot(xpts,g1_line,'r--',xpts,g2_line,'r--');
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
plot([xmin xmin axpoly(1)],[0 ymin ymin],'k--');
plot([x1 x2],[y1 y2],'ko','markersize',specs(2),'MarkerFaceColor','b');
plot(xmin,ymin,'bo',xmin,0,'bo','markersize',specs(2),'MarkerFaceColor','w'); 

subplot(nr,nc,14); hold on; 

subplot(nr,nc,15); hold on; 

subplot(nr,nc,16); hold on; axis(axpoly);
plot(xpts,g1_cube_fit,'b','linewidth',specs(1));
plot(xpts,g1_line,'r--',xpts,g2_line,'r--');
plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
plot([xmin xmin axpoly(1)],[0 ymin ymin],'k--');
plot([xmin xmin axpoly(1)],[0 chi chi],'k');
plot([x1 x2],[y1 y2],'ko','markersize',specs(2),'MarkerFaceColor','b');
plot(xmin,ymin,'bo',xmin,0,'bo','markersize',specs(2),'MarkerFaceColor','w'); 
plot(xmin,chi,'ko','markersize',specs(2),'MarkerFaceColor','b');


subplot(nr,nc,17); hold on; axis(axpoly);
plot([x1 x1 axpoly(1)],[0 ymin ymin],'k');
plot(x1,ymin,'ko','markersize',specs(2),'MarkerFaceColor','b');

%subplot(nr,nc,18); hold on; 

subplot(nr,nc,19);
semilogy(its_smooth,chifit_smooth,'r--',its(1:2),chis(1:2),'bo',...
    'linewidth',specs(1),'markersize',specs(2),'MarkerFaceColor','b');
set(gca,'xtick',xticks,'ytick',[10 100 1000 2000 6000])
xlabel(' model number (iteration)'); ylabel(' \chi(m)');
axis(ax_chi);

subplot(nr,nc,20); hold on; axis(ax_var);
plot(x_var_fit,var_red1_fit,'r--','markersize',specs(2),'linewidth',specs(1));
plot(x_var(1),var_red1(1),'bo','markersize',specs(2),'MarkerFaceColor','b');
set(gca,'xtick',xticks,'ytick',[0:20:100]);
xlabel(' model number (iteration)');
ylabel(' variance reduction');

fontsize(8); orient tall, wysiwyg

%======================================================================

break

irun0_vec = [400:20:500];
gam_vec  = [15:15:90]';
nrun = length(irun0_vec);

niter = 8;

nr=5; nc=2;
m1 = counter(1,5,1,2);
for k = 0:niter
    jj = m1(k+1);
    if mod(jj,5)==1, fontsize(7), orient tall, wysiwyg, figure; end
    disp([k jj jj*2-1 jj*2]);

    % load the chi values for this irun0 sequence
    stfm = '%3.3i';
    dir1 = '/home/carltape/wave2d/2d_adjoint/OUTPUT/';
    chi_vec = zeros(nrun,1);
    for ii=1:nrun
        irun = irun0_vec(ii) + k*2;
        stf = [dir1 'run_' num2str(sprintf(stfm,irun)) '/'];
        ww = 'summed_chi_all';
        load([stf ww '.dat']);
        chi = eval(ww);
        chi_vec(ii) = chi;
        irun_vec(ii) = irun;
    end

    subplot(nr,nc,jj*2-1); hold on;
    plot(1./gam_vec,chi_vec,'b.-','markersize',24);

    alims = axis;
    dy = alims(4) - alims(3);
    dx = alims(2) - alims(1);

    for ii=1:nrun
       irun = irun0_vec(ii) + k*2;
       text(1/gam_vec(ii),chi_vec(ii)+0.05*dy,...
           [ num2str(sprintf('%.1f', gam_vec(ii))) ' km (' num2str(irun) ')'],'fontsize',6);
       %text(1/gam_vec(ii),chi_vec(ii)+0.05*dy,[num2str(sprintf('%.1f', gam_vec(ii))) ' km'],'fontsize',6);
       %text(1/gam_vec(ii),chi_vec(ii)-0.05*dy,['(irun = ' num2str(irun_vec(ii)) ')'],'fontsize',6);
    end
    grid on;
    %xlabel(' roughness, 1 / \gamma  (\gamma determines Gaussian width for smoothing)');
    %ylabel(' misfit for m8');
    ylabel(['chi-vs-(1/\gamma) for model m' num2str(k)]);

    subplot(nr,nc,jj*2); hold on;
    plot(gam_vec,chi_vec,'b.-','markersize',24);
    for ii=1:nrun
       irun = irun0_vec(ii) + k*2;
       text(gam_vec(ii),chi_vec(ii)+0.05*dy,...
           [ num2str(sprintf('%.1f', gam_vec(ii))) ' km (' num2str(irun) ')'],'fontsize',6);
       %text(gam_vec(ii),chi_vec(ii)+0.05*dy,[num2str(sprintf('%.1f', gam_vec(ii))) ' km'],'fontsize',6);
       %text(gam_vec(ii),chi_vec(ii)-0.05*dy,['(irun = ' num2str(irun) ')'],'fontsize',6);
    end
    grid on;
    %xlabel(' smoothing, \sigma  (\sigma determines Gaussian width for smoothing)');
    %ylabel(' misfit for m8');
    ylabel(['chi-vs-\gamma for model m' num2str(k)]);
end
fontsize(7), orient tall, wysiwyg

%---------------

figure; hold on;

stc = {'r','y','g','b','m','k'};

dir1 = '/home/carltape/wave2d/2d_adjoint/gji_paper/figures/';
for ii=1:nrun
   irun0 = irun0_vec(ii);

   % fitting curve
   ww = ['chi_curve_' num2str(sprintf(stfm, irun0)) ];
   load([dir1 ww '.dat']); temp = eval(ww);
   its_smooth = temp(:,1);
   chi_smooth = temp(:,2);
   
   %ww = ['chi_' num2str(sprintf(stfm, irun0)) ];
   %load([dir1 ww '.dat']); temp = eval(ww);
   %its = temp(:,1);
   %chis = temp(:,2);
   
   plot(its_smooth,log10(chi_smooth),stc{ii});
   %plot(its,log10(chis),[stc{ii} '.'],'markersize',16);
end
grid on; xlabel(' iteration'); ylabel(' log10(chi)');
legend(['\gamma = ' num2str(gam_vec(1))],['\gamma = ' num2str(gam_vec(2))],['\gamma = ' num2str(gam_vec(3))],...
       ['\gamma = ' num2str(gam_vec(4))],['\gamma = ' num2str(gam_vec(5))],['\gamma = ' num2str(gam_vec(6))]);
for ii=1:nrun
   irun0 = irun0_vec(ii);
   ww = ['chi_' num2str(sprintf(stfm, irun0)) ];
   load([dir1 ww '.dat']); temp = eval(ww);
   its = temp(:,1);
   chis = temp(:,2);
   plot(its,log10(chis),[stc{ii} '.'],'markersize',16);
end   


if 0==1
   dir1 = '/home/store2/carltape/OUTPUT/run_001/';
   dir2 = [dir1 'event_001/'];
   
   % load source time function
   ww = 'stffor_00001_1'; load([dir2 ww]); temp = eval(ww);
   ti = temp(:,1); f = temp(:,2);
   
   % load data, synthetics, adjoint source
   ww = 'dat_00001_1'; load([dir2 ww]); temp = eval(ww); s_dat = temp(:,2);
   ww = 'syn_00001_1'; load([dir2 ww]); temp = eval(ww); s_syn = temp(:,2);
   ww = 'stfadj_00001_1'; load([dir2 ww]); temp = eval(ww); f_adj = temp(:,2);
   
   % compute velocity
   sv_dat = gradient(s_dat,ti);
   sv_syn = gradient(s_syn,ti);
   
   % specify cut time for time series
   tall = [ti f s_dat s_syn sv_dat sv_syn f_adj];
   tmax = 175;
   ikeep = find(ti <= tmax);
   tall = tall(ikeep,:);
   ti = tall(:,1);
   tall(:,7) = flipud(tall(:,7));
   
   % axes limits
   pwr_vec = [8 -4 -4 4];
   cmx_vec = [4 5 2 2];
   for ii=1:length(pwr_vec), max_vec(ii) = cmx_vec(ii) * 10^pwr_vec(ii); end
   ax_mat = [0 tmax -max_vec(1) max_vec(1) 
             0 tmax -max_vec(2) max_vec(2)
             0 tmax -max_vec(3) max_vec(3)
             0 tmax -max_vec(4) max_vec(4) ];
   
   % plot
   figure; nr=4; nc=1;
   subplot(nr,nc,1); plot(ti,tall(:,2),'b'); axis(ax_mat(1,:)); ylabel(' source time function');
   subplot(nr,nc,2); plot(ti,tall(:,4),'r--',ti,tall(:,3),'b'); axis(ax_mat(2,:)); ylabel(' displacement (data, syn)');
   subplot(nr,nc,3); plot(ti,tall(:,6),'r--',ti,tall(:,5),'b'); axis(ax_mat(3,:)); ylabel(' velocity (data, syn)');
   subplot(nr,nc,4); plot(ti,tall(:,7),'b'); axis(ax_mat(4,:)); ylabel(' tt xcorr adjoint source');
   fontsize(8); orient tall, wysiwyg
   
   % write the data to a file
   ww = ['time_series.dat'];
   fid = fopen([dir1 ww],'w');
   for ii=1:length(ti)
       fprintf(fid,'%16.7e%16.7e%16.7e%16.7e%16.7e%16.7e%16.7e\n',...
           tall(ii,1),tall(ii,2),tall(ii,3),tall(ii,4),tall(ii,5),tall(ii,6),tall(ii,7));   
   end
   fclose(fid);
   
   % write the axes info to file
   ww = ['time_series_axes.dat'];
   fid = fopen([dir1 ww],'w');
   fprintf(fid,'%16i%16i%16i%16i\n',pwr_vec(1),pwr_vec(2),pwr_vec(3),pwr_vec(4));
   fprintf(fid,'%16i%16i%16i%16i\n',cmx_vec(1),cmx_vec(2),cmx_vec(3),cmx_vec(4));
   for ii=1:length(ax_mat)
       fprintf(fid,'%16.7e%16.7e%16.7e%16.7e\n',...
           ax_mat(ii,1),ax_mat(ii,2),ax_mat(ii,3),ax_mat(ii,4));   
   end
   fclose(fid);
end

if 0==1
   dir1 = '/home/store2/carltape/OUTPUT/run_002/';   % note run 002
   dir2 = [dir1 'event_001/'];
   
   % load source time function, synthetics, adjoint source
   ww = 'stffor_00001_1'; load([dir2 ww]); temp = eval(ww); ti = temp(:,1); f = temp(:,2);
   ww = 'syn_00001_1'; load([dir2 ww]); temp = eval(ww); s_syn = temp(:,2);
   ww = 'stfadj_00001_1'; load([dir2 ww]); temp = eval(ww); f_adj = temp(:,2);
   
   % compute velocity
   sv_syn = gradient(s_syn,ti);
   
   % specify cut time for time series
   tall = zeros(length(ti),7);
   tall(:,1) = ti;
   tall(:,2) = f;
   tall(:,4) = s_syn;
   tall(:,6) = sv_syn;
   tall(:,7) = f_adj;
   tmax = 175;
   ikeep = find(ti <= tmax);
   tall = tall(ikeep,:);
   ti = tall(:,1);
   tall(:,7) = flipud(tall(:,7));
   
   % axes limits
   pwr_vec = [8 -4 -4 3];
   cmx_vec = [4 5 1.8 2];
   for ii=1:length(pwr_vec), max_vec(ii) = cmx_vec(ii) * 10^pwr_vec(ii); end
   ax_mat = [0 tmax -max_vec(1) max_vec(1) 
             0 tmax -max_vec(2) max_vec(2)
             0 tmax -max_vec(3) max_vec(3)
             0 tmax -max_vec(4) max_vec(4) ];
   
   % plot
   figure; nr=4; nc=1;
   subplot(nr,nc,1); plot(ti,tall(:,2),'b'); axis(ax_mat(1,:)); ylabel(' source time function');
   subplot(nr,nc,2); plot(ti,tall(:,4),'r--'); axis(ax_mat(2,:)); ylabel(' displacement (syn)');
   subplot(nr,nc,3); plot(ti,tall(:,6),'r--'); axis(ax_mat(3,:)); ylabel(' velocity (syn)');
   subplot(nr,nc,4); plot(ti,tall(:,7),'b'); axis(ax_mat(4,:)); ylabel(' tt xcorr adjoint source');
   fontsize(8); orient tall, wysiwyg
   
   % write the data to a file
   ww = ['time_series.dat'];
   fid = fopen([dir1 ww],'w');
   for ii=1:length(ti)
       fprintf(fid,'%16.7e%16.7e%16.7e%16.7e%16.7e%16.7e%16.7e\n',...
           tall(ii,1),tall(ii,2),tall(ii,3),tall(ii,4),tall(ii,5),tall(ii,6),tall(ii,7));   
   end
   fclose(fid);
   
   % write the axes info to file
   ww = ['time_series_axes.dat'];
   fid = fopen([dir1 ww],'w');
   fprintf(fid,'%16i%16i%16i%16i\n',pwr_vec(1),pwr_vec(2),pwr_vec(3),pwr_vec(4));
   fprintf(fid,'%16i%16i%16i%16i\n',cmx_vec(1),cmx_vec(2),cmx_vec(3),cmx_vec(4));
   for ii=1:length(ax_mat)
       fprintf(fid,'%16.7e%16.7e%16.7e%16.7e\n',...
           ax_mat(ii,1),ax_mat(ii,2),ax_mat(ii,3),ax_mat(ii,4));   
   end
   fclose(fid);
end
   
%===============================================
