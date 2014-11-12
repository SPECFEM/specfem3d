%
% function [rdist,lon,lat] = raypath(lat1,lon1,preseg,postseg,dx)
% CARL TAPE, 05-Dec-2008
% printed xxx
%
% This function returns the values of points on a great-circle "ray path",
% given a starting point and a finishing point.  The ray path is
% constructed from two segments, each extending from the starting point.
% The spacing between points will not necessarily be the same for the two
% segments, and the finishing point will not necessarily be exactly one of
% the digitized points on the ray path.
%
% calls xxx
% called by xxx
%

function [lon,lat,rdist,azis] = raypath(lon1,lat1,lon2,lat2,preseg,postseg,dx)

% distance and azimuth from point 1 to point 2
[rng, az] = distance(lat1,lon1,lat2,lon2);

% endpoints of the ray path, including the extensions
%[latA,lonA] = reckon(lat1,lon1,km2deg(preseg),wrapTo360(az+180) );
%[latB,lonB] = reckon(lat1,lon1,rng + km2deg(postseg),wrapTo360(az+180) );

% lengths of two segments and number of legs per segment
dist1 = preseg;
dist2 = postseg+deg2km(rng);
n1 = ceil(dist1/dx)+1;
n2 = ceil(dist2/dx)+1;

% two segments
[lats1,lons1] = track1(lat1,lon1,wrapTo360(az+180),km2deg(dist1),[],[],n1);
[lats2,lons2] = track1(lat1,lon1,az,km2deg(dist2),[],[],n2);
npt1 = length(lats1);
npt2 = length(lats2);
npt = npt1 + npt2;

% connect everything
lon = zeros(npt,1);
lat = zeros(npt,1);
lon = [flipud(lons1) ; lons2(2:end) ];
lat = [flipud(lats1) ; lats2(2:end) ];
n = length(lon);

[ddist,azis] = distance(lat1*ones(n,1),lon1*ones(n,1),lat,lon);
rdist = deg2km(ddist);
rdist(1:npt1) = -rdist(1:npt1);

%------------------------------------

% EXAMPLE
if 0==1
    lon1 = -120; lat1 = 36;
    lon2 = -117; lat2 = 33;
    preseg = 50; postseg = 50;
    dx = 10;
    
    [lon,lat,rdist,azis] = raypath(lon1,lat1,lon2,lat2,preseg,postseg,dx);
    figure; hold on;
    plot(lon,lat,'.');
    plot(lon1,lat1,'rp',lon2,lat2,'ko','markersize',16);
    for ii=1:length(lon), text(lon(ii),lat(ii),num2str(ii)); end
    
    for ii=1:length(lon)
       disp(sprintf('lon, lat, dist, az: %10.3f %10.3f %10.3f %10.3f',...
           lon(ii),lat(ii),rdist(ii),azis(ii)));
    end
end

%======================================================================
