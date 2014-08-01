function plot_cc(data,syn,new,new_dt,dlnA_cc,dt_cc)

fsize = 11;
lsize = 1.0;
msize = 12;
cgreen = [0.1210 0.6050 0];
f1 = 0; f2 = 0.25;

% strings for cross-correlation measurements
stdT_cc = sprintf(' dT-cc = %.3f +/- %.3f',dt_cc(1),dt_cc(2));
stdA_cc = sprintf(' dlnA-cc = %.3f +/- %.3f',dlnA_cc(1),dlnA_cc(2));

subplot(2,1,1), hold on
%title('Am','fontsize',fsize)
plot(data(:,1),data(:,2),'-k','linewidth',lsize)
plot(syn(:,1),syn(:,2),'--r','linewidth',lsize)
plot(new_dt(:,1),new_dt(:,2),'--b','linewidth',lsize)
plot(new(:,1),new(:,2),'-','linewidth',lsize,'Color',cgreen)
xlabel('Time [s]','fontsize',fsize)
% HOW CAN WE GET THE LEGEND TO PLOT WHERE WE WANT IT TO?
% MY QUICK FIX WAS TO PUT IT OUTSIDE, BUT I HAD TO SHORTEN THE LABELS.
legend('Data','Syn-0','Syn-1','Syn-2','location','NorthWestOutside')

subplot(2,2,3), hold on
plot_bars(f1,f2,dlnA_cc(1),dlnA_cc(2),cgreen,lsize)     % cc measurement
ylabel('\Delta lnA','fontsize',fsize)
xlabel('Frequency [Hz]','fontsize',fsize)
  title({'Amplitude measurement',stdA_cc},'fontsize',fsize)

subplot(2,2,4), hold on
plot_bars(f1,f2,dt_cc(1),dt_cc(2),cgreen,lsize)     % cc measurement
xlabel('Frequency [Hz]','fontsize',fsize)
ylabel('\Delta T [s]','fontsize',fsize)
title({'Traveltime measurement',stdT_cc},'fontsize',fsize)

