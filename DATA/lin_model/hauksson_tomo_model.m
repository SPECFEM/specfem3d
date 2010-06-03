%
% hauksson_tomo_model.m
% Carl Tape, 11-July-2007
% 
% This program reads in the southern California tomography models of
%    Hauksson (JGR 2000)
%    Lin-Shearer-Hauksson-Thurber (JGR 2007)
% and outputs some figures.
%
% NOTE: THIS PROGRAM REQUIRES SOME ADDITIONAL CUSTOM SUBROUTINES.
%       If you want to run the program and generate the output
%       files and figures, email Carl Tape for the additional scripts.
%
% calls xxx
% called by xxx
%

clear
close all
clc
format compact

idir = '/home/carltape/gmt/tomography/';

ihauk = 1;    % =1 to also plot Hauksson (2000) model
iwrite = 0;

% bounds for 'standard' SPECFEM3D simulation
% LATITUDE_MIN                    = 32.2d0
% LATITUDE_MAX                    = 36.8d0
% LONGITUDE_MIN                   = -120.3d0
% LONGITUDE_MAX                   = -114.7d0
axbox = [-120.3 -114.7 32.2 36.8];

% bounds for Hauksson regular grid
%   double precision, parameter :: ORIG_LONG_HAUKSSON = -121.d0,END_LONG_HAUKSSON = -114.d0
%   double precision, parameter :: ORIG_LAT_HAUKSSON = 32.d0,END_LAT_HAUKSSON = 37.d0
axbox_hauk = [-121 -114 32 37];

%------------------------------------------------------
% dimensions of PoChen et al. LA Basin model (2007)
%
% ---------- Forwarded message ----------
% Date: Tue, 30 Oct 2007 14:33:52 -0400
% From: Po Chen <pseudopochen@gmail.com>
% To: Carl Tape <carltape@gps.caltech.edu>
% Subject: Re: BSSA paper, figure 1 bounds
% 
% Hello Carl,
% 
% I could not find the exact lat/lon for that box used in the paper. The
% following numbers might be close:
% 
% lon: -118.76 -117.22
% lat: 33.66 34.41
% 
% Yes, the depth extend of the inversion is 26km.
% 
% Best,
% Po
if 0==1
    ax0 = [-118.76 -117.22 33.66 34.41];
    ax_utm = axes_utm2ll(ax0,11,0);
    dx = diff(ax_utm(1:2));
    dy = diff(ax_utm(3:4));
    dz = 26*1e3;
    V = dx*dy*dz;
    disp(sprintf(' dx dy dz : %.2f km %.2f km %.2f km',dx*1e-3,dy*1e-3,dz*1e-3));
    disp(sprintf(' Volume : %.4e m^3',V));
    
    lons = [238.6100 239.0480 239.4235 238.9854];
    lats = [35.6421 35.2860 35.5913 35.9474];
    disp([lons'-360 lats']);
    
    % bounds of our simulations
    ax0 = [-121.6 -114.7 32.2 36.8];
    ax0_utm = axes_utm2ll(ax0,11,0);
    dx0 = diff(ax0_utm(1:2));
    dy0 = diff(ax0_utm(3:4));
    dz0 = 60*1e3;
    V0 = dx0*dy0*dz0;
    disp(sprintf(' dx dy dz : %.2f km %.2f km %.2f km',dx0*1e-3,dy0*1e-3,dz0*1e-3));
    disp(sprintf(' Volume : %.4e m^3',V0));
    disp(sprintf(' Volume ratio : %.2f m^3',V0/V));
    break
end
%------------------------------------------------------

if 0==1
    % regrid_hauksson_regular.f90
    % write(13,'(2i10,3e20.10)') i, j, utm_x_new(i,j), utm_y_new(i,j), distmin_all(i,j)
    dall = load('distmin.dat');
    
    utmx = dall(:,3);
    utmy = dall(:,4);
    dmin = dall(:,5)/1000;    % km
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on;
    plot3(utmx,utmy,zeros,'r.','markersize',1);
    plot3(utmx,utmy,dmin,'b.','markersize',1); grid on; box on;
    xlabel('UTM-X (m)'); ylabel('UTM-Y (m)'); zlabel(' Min distance, km');
    subplot(nr,nc,2); plot_histo(dmin,[0:1:10]); xlabel(' Min distance, km');
    orient tall, wysiwyg
    
    break
end

%------------------------------------------------------
% Hauksson (2000) tomo model

if ihauk==1
    disp('  '); disp(' Hauksson (2000) tomography model');

    dir_hauk = [idir 'hauk_2000/'];
    lines = textread([dir_hauk 'hauksson_new_format.dat'],'%s','delimiter','\n','whitespace','');
    nlines = length(lines);

    % read the file into a vector
    vec_all = [];
    for ii = 1:nlines
        a = str2num(lines{ii});
        vec_all = [vec_all ; a(:)];
    end
    nlen = length(vec_all);
    clear lines

    % depth levels (constants.h)
    udep_hauk = [1 4 6 10 15 17 22 31 33]';
    ndep_hauk = length(udep_hauk);

    % number of entries per lat-lon point
    npts_per_latlon = 2 + 2*ndep_hauk;
    npts_hauk = nlen / npts_per_latlon;   % 27 * 41 = 1107
    NX_hauk = 27;
    NY_hauk = 41;

    % convert this into a format similar to Lin2007
    % LON, LAT, DEP, Vp, Vs
    %npts_hauk = 2;
    hauk_model = zeros(npts_hauk*ndep_hauk,5);
    kk = 1;
    jj = 1;
    for ii = 1:npts_hauk
        % extract a 'block' of values
        kmin = kk;
        kmax = kk + npts_per_latlon - 1;
        temp = vec_all(kmin:kmax);

        % fill in a new matrix
        jmin = jj;
        jmax = jmin + ndep_hauk - 1;
        hauk_model(jmin:jmax,:) = ...
          [ temp(1)*ones(ndep_hauk,1) temp(2)*ones(ndep_hauk,1) ...
            udep_hauk  temp(3:11)  temp(12:20) ];

        kk = kmax + 1;
        jj = jmax + 1;
    end

    lon_hauk = hauk_model(:,1);
    lat_hauk = hauk_model(:,2);
    dep_hauk = hauk_model(:,3);
    
    % find the corners (clockwise from the top)
    [i1, i2] = max(lat_hauk); lon1 = lon_hauk(i2); lat1 = lat_hauk(i2);
    [i1, i2] = max(lon_hauk); lon2 = lon_hauk(i2); lat2 = lat_hauk(i2);
    [i1, i2] = min(lat_hauk); lon3 = lon_hauk(i2); lat3 = lat_hauk(i2);
    [i1, i2] = min(lon_hauk); lon4 = lon_hauk(i2); lat4 = lat_hauk(i2);
    disp(' Corners of Hauksson (2000) model : ');
    st1ll = sprintf('(%8.2f,%7.2f)',lon1,lat1); disp(st1ll);
    st2ll = sprintf('(%8.2f,%7.2f)',lon2,lat2); disp(st2ll);
    st3ll = sprintf('(%8.2f,%7.2f)',lon3,lat3); disp(st3ll);
    st4ll = sprintf('(%8.2f,%7.2f)',lon4,lat4); disp(st4ll);
    disp(['Grid is ' num2str(NX_hauk) ' by ' num2str(NY_hauk) ' = ' num2str(npts_hauk) ' nodes']);
    disp([ num2str(ndep_hauk) ' depth layers']);

    % extract the upper layer
    isurface_hauk = find(dep_hauk == udep_hauk(1));
    
    %-----------------------------------------------
    % load final smoothed version of Hauksson model

end

%------------------------------------------------------
% Lin-Shearer-Hauksson-Thurber (2007) tomography model

disp('  '); disp(' Lin-Shearer-Hauksson-Thurber (2007) tomography model');

% load Vp model
dir_lin = [idir 'lin_2007/'];
vp = load([dir_lin 'pa_all-1']);
alpha = vp(:,1);
lon = vp(:,3);
lat = vp(:,2);
gridX = vp(:,4);
gridY = vp(:,5);
dep = vp(:,6);
n = length(alpha);

% load Vp/Vs model
vpvs = load([dir_lin 'ra_all-1']);
beta = 1./vpvs(:,1) .* alpha;

% find the depth levels
udep = unique(dep);
ndep = length(udep);

% assuming a max depth of the model, compute the thickness of each layer
dmax = 60;
dlayer = diff(udep)/2;
dtop = udep - [0; dlayer];
dbot = udep + [dlayer ; dmax-udep(end)];
dthickness = dbot-dtop;
dall = [udep dtop dbot dthickness];
if sum(dall(:,4)) ~= dmax, error(['thickness of all layers should equal ' num2str(dmax) ]); end
disp('Here is how we might think of the depth layers:');
disp('   reference    top    bottom   thickness');
disp(dall);

% find range of grid
mingY = min(gridY); maxgY = max(gridY);
mingX = min(gridX); maxgX = max(gridX);
rnY = maxgY - mingY;
rnX = maxgX - mingX;

% find the corners (clockwise from the top)
[i1, i2] = max(lat); y1 = gridY(i2); x1 = gridX(i2); lon1 = lon(i2); lat1 = lat(i2);
[i1, i2] = max(lon); y2 = gridY(i2); x2 = gridX(i2); lon2 = lon(i2); lat2 = lat(i2);
[i1, i2] = min(lat); y3 = gridY(i2); x3 = gridX(i2); lon3 = lon(i2); lat3 = lat(i2);
[i1, i2] = min(lon); y4 = gridY(i2); x4 = gridX(i2); lon4 = lon(i2); lat4 = lat(i2);
st1ll = sprintf('(%8.2f,%7.2f)',lon1,lat1);
st2ll = sprintf('(%8.2f,%7.2f)',lon2,lat2);
st3ll = sprintf('(%8.2f,%7.2f)',lon3,lat3);
st4ll = sprintf('(%8.2f,%7.2f)',lon4,lat4);
st1 = sprintf('(%5i,%5i)',x1,y1);
st2 = sprintf('(%5i,%5i)',x2,y2);
st3 = sprintf('(%5i,%5i)',x3,y3);
st4 = sprintf('(%5i,%5i)',x4,y4);
disp(' Corners of Lin-Shearer-Hauksson-Thurber (2007) model : ');
disp(st1ll); disp(st2ll); disp(st3ll); disp(st4ll);
disp(st1); disp(st2); disp(st3); disp(st4);

% extract the upper layer
isurface = find(dep == udep(1));
izeroY = find(and(gridX == 0, dep == udep(1)));
izeroX = find(and(gridY == 0, dep == udep(1)));
iorigin = find(and(and(gridY == 0, gridX == 0), dep == udep(1)) );
NX = length(izeroX);
NY = length(izeroY);
npts = NX * NY;
disp(['Grid is ' num2str(NX) ' by ' num2str(NY) ' = ' num2str(npts) ' nodes']);
disp([ num2str(ndep) ' depth layers']);

% extract the boundary lat-lon points (for GMT plots)
% ORDERING: from top point, clockwise
iNE_boundary = flipud( find( gridX(isurface) == mingX ) );
iSE_boundary =         find( gridY(isurface) == mingY );
iSW_boundary =         find( gridX(isurface) == maxgX );
iNW_boundary = flipud( find( gridY(isurface) == maxgY ) );
iboundary = isurface([iNE_boundary ; iSE_boundary(2:end) ; iSW_boundary(2:end) ; iNW_boundary]);
%figure; plot(lon(iboundary),lat(iboundary),'b.-');
if iwrite == 1
    write_xy_points([dir_lin 'lin_boundary'],lon(iboundary),lat(iboundary));
end

% compute the mean 1D model -- including scaling to density from alpha
alpha_1D = zeros(ndep,3);
beta_1D  = zeros(ndep,3);
for ii=1:ndep
    inds = find(dep == udep(ii));
    alpha_1D(ii,:) = [min(alpha(inds))  mean(alpha(inds)) max(alpha(inds))] ;
    beta_1D(ii,:)  = [min(beta(inds))   mean(beta(inds))  max(beta(inds)) ] ;
end
rho_1D = alpha_rho(alpha_1D*1e3);
disp('  ');
disp('1D averaged model:');
disp('     depth  thickness   vp-min   vp-mean    vp-max    vs-min   vs-mean    vs-max    rho-min  rho-mean  rho-max');
disp([ udep dthickness alpha_1D beta_1D rho_1D/1000]);

% using the 1D model and the layer thicknesses, compute the overall mean velocities
alpha_mean = sum(alpha_1D(:,2) .* dthickness) / sum(dthickness);
beta_mean  = sum(beta_1D(:,2) .* dthickness) / sum(dthickness);
c_mean = sqrt( alpha_mean^2 - (4/3)*beta_mean^2 );
disp('  ');
disp([' Overall mean velocities, assuming a bottom depth of ' num2str(dmax) ' km:'])
disp(sprintf('%8s : %.1f (%.4f) km/s','alpha',alpha_mean,alpha_mean));
disp(sprintf('%8s : %.1f (%.4f) km/s','beta',beta_mean,beta_mean));
disp(sprintf('%8s : %.1f (%.4f) km/s','c',c_mean,c_mean));

%---------------------------------
% figures

[xmat, ymat] = horzlines(-udep,0,1);

figure; hold on;
plot(xmat,ymat,'b','linewidth',2);
if ihauk==1
    [xmat, ymat] = horzlines(-udep_hauk,0,1);
    plot(xmat,ymat,'r--','linewidth',2);
    title({['BLUE : Lin et al. 2007 (' num2str(ndep) ' layers)'],
        ['RED : Hauksson 2000 (' num2str(ndep_hauk) ' layers)']});
end
ylim([-35 5]); ylabel(' Depth (km)');

figure; hold on;
plot(lon(isurface),lat(isurface),'b.');
if ihauk==1
plot(lon_hauk(isurface_hauk),lat_hauk(isurface_hauk),'r.');
end
axis equal; ax0 = axis; axis(axes_expand(ax0,1.05,1));
plot(axbox_hauk([1 2 2 1 1]),axbox_hauk([3 3 4 4 3]),'r','linewidth',2);
plot(axbox([1 2 2 1 1]),axbox([3 3 4 4 3]),'k','linewidth',2);
if ihauk==1, legend('Lin et al. 2007','Hauksson 2000','Regular tomo in SPECFEM','SPECFEM bounds');
else legend('Lin et al. 2007','Regular tomo in SPECFEM','SPECFEM bounds'); end
xlabel(' Longitude'); ylabel(' Latitude');
title(' Nodes for SoCal tomography model');
orient tall, wysiwyg

figure; hold on;
plot(lon(isurface),lat(isurface),'k.')
axis equal; ax0 = axis; axis(axes_expand(ax0,1.05,1));
plot(lon(izeroY),lat(izeroY),'ro');
plot(lon(izeroX),lat(izeroX),'bo');
text(lon1,lat1,st1); text(lon2,lat2,st2); text(lon3,lat3,st3); text(lon4,lat4,st4);
plot(lon(iorigin),lat(iorigin),'kp','markersize',15);
xlabel(' Longitude'); ylabel(' Latitude');
title(' Nodes for Lin et al. (2007) model');
orient tall, wysiwyg

figure; hold on;
plot(gridX(isurface),gridY(isurface),'k.')
axis equal; ax0 = axis; axis(axes_expand(ax0,1.05,1));
plot(gridX(izeroY),gridY(izeroY),'ro');
plot(gridX(izeroX),gridY(izeroX),'bo');
plot(gridX(iorigin),gridY(iorigin),'kp','markersize',15);
text(x1,y1,st1ll); text(x2,y2,st2ll); text(x3,y3,st3ll); text(x4,y4,st4ll);
xlabel(' NX, grid coordinates (NE -- SW)  (km)');
ylabel(' NY, grid coordinates (SE -- NW)  (km)');
title({' Nodes for Lin et al. (2007) model',
    'Looking from BELOW -- NORTH points to the upper left'});
orient tall, wysiwyg

figure; hold on;
plot3(lon,lat,-dep,'b.');
if ihauk==1
plot3(lon_hauk,lat_hauk,-dep_hauk,'r.');
legend('Lin et al. 2007','Hauksson 2000');
end
view(3); box on, grid on;
xlabel(' Longitude'); ylabel(' Latitude'); zlabel(' Depth (km)');
orient tall, wysiwyg

if iwrite == 1
    file = 'lin_new_format.dat';
    fid = fopen([dir_lin file],'w');
    inds0 = [1 : npts : n ];
    for kk = 1:npts
        inds = inds0 + (kk-1);
        atemp = alpha(inds);
        btemp =  beta(inds);
        
        fprintf(fid,'%12.6f%12.6f%12.6f\n',lon(kk),lat(kk),atemp(1));
        fprintf(fid,'%12.6f%12.6f%12.6f\n',atemp(2),atemp(3),atemp(4));
        fprintf(fid,'%12.6f%12.6f%12.6f\n',atemp(5),atemp(6),atemp(7));
        fprintf(fid,'%12.6f%12.6f%12.6f\n',atemp(8),btemp(1),btemp(2));
        fprintf(fid,'%12.6f%12.6f%12.6f\n',btemp(3),btemp(4),btemp(5));
        fprintf(fid,'%12.6f%12.6f%12.6f\n',btemp(6),btemp(7),btemp(8));
    end
    fclose(fid);
end

%==========================================
