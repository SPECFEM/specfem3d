clear 

% $$$ MESHER prem_ani
% $$$ period | max  |
% $$$  T[s]  |nproc | memory [GB]
% $$$ -------|------|-------------
% $$$ 
% $$$  50.0  |   8  |    < 1
% $$$  25.0  |  12  |    < 1
% $$$  20.0  |  16  |    < 1
% $$$  15.0  |  20  |    < 1
% $$$  10.0  |  28  |    < 1
% $$$   9.0  |  32  |    < 1
% $$$   5.0  |  56  |    < 1
% $$$   4.4  |  64  |    < 1
% $$$   4.0  |  68  |    < 1
% $$$   3.0  |  92  |    2.0
% $$$   2.2  | 128  |    3.5
% $$$   2.0  | 136  |    4.0
% $$$   1.1  | 256  |   13.0
% $$$   1.0  | 272  |   16.0

mesher_period = [50 25 20 15 10 9 5 4.4 4 3 2.2 2 1.1 1];
mesher_ram = [1 1 1 1 1 1 1 1 1 2 3.5 4 13 16];
mesher_maxproc = [8 12 16 20 28 32 56 64 68 92 128 136 256 272];

figure(5);clf
[ax,h1,h2]=plotyy(mesher_period,mesher_ram,mesher_period,mesher_maxproc,'loglog')
%hold on
set(h2,'LineWidth',2)
set(h1,'LineWidth',2)
set(gca,'FontSize',20,'XTick',[1 10 20 50])
grid on
axtmp = findobj(gcf,'Type','axes');
axes(axtmp(2))
ylabel('RAM memory [GB]')
set(gca,'FontSize',20)
axes(axtmp(1))
set(gca,'FontSize',20)
ylabel('max # cores')
set(gca,'FontSize',20)
legend('cores','memory')
set(gca,'FontSize',20)
title('Mesher: computational constraints')
xlabel('seismic period [s]')
axis([1 70 1 300])
set(axtmp(2),'FontSize',20)
set(axtmp(2),'YTick',[1 10 20])
axes(axtmp(1))
set(axtmp(1),'YTick',[16 64 256])

% $$$ 
% $$$ SOLVER prem_ani
% $$$ total cpu time (time * nproc) in time loop for 1000s simulation
% $$$ coarse grained attenuation with 5 linear solids
% $$$ dt as suggested by the mesher
% $$$ 
% $$$ Monte Rosa: cray Xe6 with AMD Interlagos, gfortran compiler:
% $$$ 
% $$$ period |       |       |    total cpu time - [s]
% $$$  T[s]  | dt[s] | nproc | monopole | dipole | quadpole
% $$$ -------|-------|-------|----------|--------|---------
% $$$  50.0  | 0.12  |   4   |     280  |   390  |   450
% $$$  25.0  | 0.12  |   4   |    1000  |  1360  |  1460
% $$$  15.0  | 0.12  |   8   |    2930  |  3920  |  4110
% $$$  10.0  | 0.11  |   8   |    7120  |  8820  |  9300
% $$$   5.0  | 0.056 |  16   |   57700  |        |
% $$$   3.0  | 0.033 |  32   |  257072  |        |
% $$$ 

solver_period = [50 25 15 10 5 3];
dt = [0.12 0.12 0.12 0.11 0.056 0.033];
solver_proc = [4 4 8 8 16 32];
cpu_mono=[280 1000 2930 7120 57700 257072];
cpu_mono_30min_hrs=cpu_mono*1.8/60./60.;

figure(6);clf
h=loglog(solver_period,cpu_mono_30min_hrs)
set(h,'LineWidth',2)
set(gca,'FontSize',20,'XTick',[1 2 5 10 20 50],'YTick',[0.1 1 10 100])
grid on
ylabel('Total CPU time [hrs]')
set(gca,'FontSize',20)
title('30min seismogram: Total CPU cost')
xlabel('seismic period [s]')
axis([1 70 0 150])

figure(7);clf
h=loglog(solver_period,cpu_mono_30min_hrs./solver_proc,'x-')
hold on
for i=1:length(solver_proc)
h1=text(solver_period(i)*1.1,cpu_mono_30min_hrs(i)/solver_proc(i),sprintf('%i', ...
                                                  solver_proc(i)))
set(h1,'FontSize',20)
end
set(h,'LineWidth',2,'MarkerSize',10)
set(gca,'FontSize',20,'XTick',[1 2 5 10 20 50],'YTick',[0.01 0.1 1 5])
grid on
ylabel('Wall-clock CPU time [hrs]')
set(gca,'FontSize',20)
title('30min seismogram: parallel job length')
xlabel('seismic period [s]')
%axis([1 70 1 10])

% $$$ 
% $$$ my office computer (Intel Quad Core) is about 3 times faster (per core)
% $$$ 
% $$$ period |       |
% $$$  T[s]  | dt[s] | monopole
% $$$ -------|-------|----------
% $$$  25.0  | 0.12  |    331  |
% $$$  15.0  | 0.12  |    866  |
% $$$  10.0  | 0.11  |   2100  |
% $$$   5.0  | 0.056 |  28000  |
  
% $$$ STRONG SCALING
% $$$  - prem light, 2.1s
% $$$  - auf 8, 16, 32, 64 und 128 cores
% $$$  - 1000 zeitschritte
% $$$  - anzahl elemente nimmt leicht zu bei mehr cores:
% $$$ 
% $$$ #cores
% $$$     timeloop time in s
% $$$  		 nel per proc
% $$$   8 2569.837158  254760
% $$$  16 1122.863037  125732
% $$$  32  535.962036   63792
% $$$  64  263.186005   32500
% $$$ 128  149.062012   18396
% $$$ 

% SCALING
nproc_strong = [8 16 32 64 128];
time=[2569.837158 1122.863037 535.962036 263.186005  149.062012];
numel=[254760 125732 63792 32500 18396];
time_theor=2569.837158*8./nproc_strong;

figure(1);clf
h1=loglog(nproc_strong,time,'r-x')
grid on
set(h1,'LineWidth',2,'MarkerSize',14)
hold on
h2=loglog(nproc_strong,time_theor,'k--')
set(gca,'FontSize',20,'YTick',[200 500 1000 2000],'XTick',[8 16 32 ...
                    64 128])
legend('AxiSEM','optimal')
title('Strong scaling')
xlabel('# cores')
ylabel('CPU time for 1000 time steps [s] ]')
set(h2,'LineWidth',2)
axis([0 150 0 3000])

% $$$ WEAK SCALING
% $$$ #nproc
% $$$      #runtime	 #elems/proc
% $$$ 004  240.111008  35428
% $$$ 008  319.252014  35664
% $$$ 016  302.423004  35568
% $$$ 032  302.613007  35688
% $$$ 064  302.783020  35600
% $$$ 128  314.622009  35724
% $$$ 256  329.309021  36300
nproc_weak = [8 16 32 64 128 256];
time_weak=[319.252014 302.423004 302.613007 302.783020 314.622009 329.309021];
numel_weak=[35664 35568 35688 35600 35724 36300];
time_theor_weak=302.423004.*nproc_weak./nproc_weak;

figure(2);clf
h1=semilogx(nproc_weak,time_weak,'r-x')
set(h1,'LineWidth',2,'MarkerSize',14)
hold on
grid on
h2=semilogx(nproc_weak,time_theor_weak,'k--')
set(gca,'FontSize',20,'XTick',[8 16 32 ...
                    64 128 256])
legend('AxiSEM','optimal')
title('Weak scaling')
xlabel('# cores')
ylabel('CPU time for 1000 time steps [s] ]')
set(h2,'LineWidth',2)
axis([0 300 100 500 ])