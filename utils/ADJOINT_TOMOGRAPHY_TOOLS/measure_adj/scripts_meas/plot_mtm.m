function plot_mtm(data,syn,dataf,synf,new,new_dt,dlnA,dt,dlnA_full,dt_full,dlnA_ave,dt_ave,...
		  new_cc,new_cc_dt,dlnA_cc,dt_cc,flims,iall)

fsize = 9;
lsize = 1.0;
msize = 12;

cgreen = [0.1210 0.6050 0];

fmax = max(flims);
f1 = -0.05*fmax;
f2 = 1.05*fmax;

% strings for cross-correlation measurements
stdT_cc = sprintf(' dT-cc (green) = %.3f +/- %.3f',dt_cc(1),dt_cc(2));
stdA_cc = sprintf(' dlnA-cc (green) = %.3f +/- %.3f',dlnA_cc(1),dlnA_cc(2));
stdT_mtave = sprintf(' dT-mt-ave (blue) = %.3f +/- %.3f',dt_ave(1),dt_ave(2));
stdA_mtave = sprintf(' dlnA-mt-ave (blue) = %.3f +/- %.3f',dlnA_ave(1),dlnA_ave(2));

subplot(3,1,1), hold on
%title('Am','fontsize',fsize)
plot(data(:,1),data(:,2),'-','linewidth',2*lsize,'Color','k')
plot(syn(:,1),syn(:,2),'-','linewidth',2*lsize,'Color','r')

if 0==1
  plot(new_dt(:,1),new_dt(:,2),'--','linewidth',lsize,'Color','cyan')
  plot(new(:,1),new(:,2),'--','linewidth',lsize,'Color','cyan')
  plot(new_cc_dt(:,1),new_cc_dt(:,2),'--','linewidth',lsize,'Color',cgreen)
  plot(new_cc(:,1),new_cc(:,2),'--','linewidth',lsize,'Color',cgreen)
  legend('Data','Syn','Syn-MT-1','Syn-MT-2','Syn-CC-1','Syn-CC-2','location','BestOutside')
else
  plot(new(:,1),new(:,2),'--','linewidth',lsize,'Color','cyan')
  plot(new_cc(:,1),new_cc(:,2),'--','linewidth',lsize,'Color',cgreen)
  legend('Data','Syn','Syn-MT','Syn-CC','location','BestOutside')
end
xlabel('Time [s]','fontsize',fsize)

% plot power of data and synthetics
subplot(3,3,4), hold on
plot(dataf(:,1),dataf(:,2),'-','linewidth',2*lsize,'Color','k')
plot(synf(:,1),synf(:,2),'-','linewidth',2*lsize,'Color','r')
xlim([f1 f2])
ax2 = axis; [xmat, ymat] = vertlines(flims,ax2(3),ax2(4)); plot(xmat(:,1:4),ymat(:,1:4),'k',xmat(:,5:6),ymat(:,5:6),'r');
ylabel('Power','fontsize',fsize)
xlabel('Frequency [Hz]','fontsize',fsize)
title('Power of data and synthetics','fontsize',fsize)

subplot(3,3,5), hold on
% plot multitaper measurements
if iall == 1
errorbar(dlnA(:,1),dlnA(:,2),dlnA(:,3),'.','markersize',msize,'linewidth',lsize)
ax1 = axis;
plot(dlnA_full(:,1),dlnA_full(:,2),'b','linewidth',lsize/2)
plot(dlnA_full(:,1),dlnA_full(:,2)+dlnA_full(:,3),'r--','linewidth',lsize/2)  % top error
plot(dlnA_full(:,1),dlnA_full(:,2)-dlnA_full(:,3),'r--','linewidth',lsize/2)  % bottom error
axis([f1 f2 ax1(3:4)])
end
plot_bars(f1,f2,dlnA_ave(1),dlnA_ave(2),[0 0 1],lsize)  % mt average measurement
plot_bars(f1,f2,dlnA_cc(1),dlnA_cc(2),cgreen,lsize)     % cc measurement
ax2 = axis; [xmat, ymat] = vertlines(flims,ax2(3),ax2(4)); plot(xmat(:,1:4),ymat(:,1:4),'k',xmat(:,5:6),ymat(:,5:6),'r');
ylabel('\Delta lnA','fontsize',fsize)
xlabel('Frequency [Hz]','fontsize',fsize)
  title({'Amplitude measurement',stdA_cc,stdA_mtave},'fontsize',fsize)

subplot(3,3,6), hold on
% plot multitaper measurements
if iall==1
errorbar(dt(:,1),dt(:,2),dt(:,3),'.','markersize',msize,'linewidth',lsize)
ax1 = axis;
plot(dt_full(:,1),dt_full(:,2),'b','linewidth',lsize/2)
plot(dt_full(:,1),dt_full(:,2)+dt_full(:,3),'r--','linewidth',lsize/2)  % top error
plot(dt_full(:,1),dt_full(:,2)-dt_full(:,3),'r--','linewidth',lsize/2)  % bottom error
axis([f1 f2 ax1(3:4)])
end
plot_bars(f1,f2,dt_ave(1),dt_ave(2),[0 0 1],lsize)  % mt average measurement
plot_bars(f1,f2,dt_cc(1),dt_cc(2),cgreen,lsize)     % cc measurement
ax2 = axis; [xmat, ymat] = vertlines(flims,ax2(3),ax2(4)); plot(xmat(:,1:4),ymat(:,1:4),'k',xmat(:,5:6),ymat(:,5:6),'r');
xlabel('Frequency [Hz]','fontsize',fsize)
ylabel('\Delta T [s]','fontsize',fsize)
title({'Traveltime measurement',stdT_cc,stdT_mtave},'fontsize',fsize)

