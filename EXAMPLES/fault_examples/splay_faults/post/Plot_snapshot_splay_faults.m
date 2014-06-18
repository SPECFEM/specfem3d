addpath ~/SPECFEM3D/Post-processing

clear all 
close all
ITSNAP =   33000;
XLIM   = [0 120]; % range along-strike
SLIM   =  [-8 8];    % slip range

dA  = FSEM3D_snapshot(ITSNAP,'.',1);
dBC = FSEM3D_snapshot(ITSNAP,'.',2);

% along-dip slip is negative for thrust faulting
dA.Dz =   -dA.Dz;
dBC.Dz = -dBC.Dz;

% shift to Wendt's coordinate system
dA.Y = 180-dA.Y;
dA.X = 227-dA.X;
dBC.Y = 180-dBC.Y;
dBC.X = 90-dBC.X;

figure(1)
subplot(221)
plot(dA.Y,dA.Dx,'.')
set(gca,'XLim',XLIM)
ylabel('Dx (m)')
subplot(223)
plot(dA.Y,dA.Dz,'.')
set(gca,'XLim',XLIM)
ylabel('Dz (m)')
xlabel('Y along-strike (km)')
subplot(222)
plot(dBC.Y,dBC.Dx,'.')
set(gca,'XLim',XLIM)
subplot(224)
plot(dBC.Y,dBC.Dz,'.')
set(gca,'XLim',XLIM)

figure(2)
subplot(121)
[hc,hf]=plotclr(dA.Y,dA.X,dA.Dz,'s',SLIM);
set(hf,'YDir','reverse')
xlabel('Along strike (km)')
ylabel('Down dip (km)')
axis equal
axis tight
set(hf,'XLim',XLIM)

subplot(122)
[hc,hf]=plotclr(dBC.Y,dBC.X,dBC.Dz,'s',SLIM);
set(hf,'YDir','reverse')
xlabel('Along strike (km)')
ylabel('Down dip (km)')
axis equal
axis tight
set(hf,'XLim',XLIM)

% This is the colormap of Wendt et al, but it does not work well yet with plotclr:
%colormap('jet');
%map1 = colormap;
%l1 = length(map1);
%map2 = map1(l1/2+1:end,:);
%colormap(map2);

