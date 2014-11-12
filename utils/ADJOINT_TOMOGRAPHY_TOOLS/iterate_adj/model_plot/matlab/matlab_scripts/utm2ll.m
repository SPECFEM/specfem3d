%
% function [x,y,estr] = utm2ll(xi,yi,s_zone,i_type) 
%
% Black-box Matlab code to convert between utm and lon-lat coordinates.
% NOTE: Requires Matlab Mapping Toolbox.
%
% To find a UTM zone for a specific point, use utmzone(lat,lon)
% UTM zones: http://www.dmap.co.uk/utmworld.htm
%
% Ellipsoids available: "help almanac"
% everest       1830 Everest ellipsoid
% bessel        1841 Bessel ellipsoid
% airy          1849 Airy ellipsoid
% clarke66      1866 Clarke ellipsoid
% clarke80      1880 Clarke ellipsoid
% international 1924 International ellipsoid
% krasovsky     1940 Krasovsky ellipsoid
% wgs60         1960 World Geodetic System ellipsoid
% iau65         1965 International Astronomical Union ellipsoid
% wgs66         1966 World Geodetic System ellipsoid
% iau68         1968 International Astronomical Union ellipsoid
% wgs72         1972 World Geodetic System ellipsoid
% grs80         1980 Geodetic Reference System ellipsoid
% wgs84         1984 World Geodetic System ellipsoid
%
% Mapping Toolbox ellipsoid representations are two-element vectors,
% called ellipsoid vectors. The ellipsoid vector has the form
% [semimajor_axis eccentricity]. The semimajor axis can be in any
% unit of distance; the choice of units typically drives the units
% used for distance outputs in the toolbox functions.
% Meters, kilometers, or Earth radii (i.e., a unit sphere) are most
% frequently used. See Functions that Define Ellipsoid Vectors for details.
%
% calls xxx
% called by test_utm.m
%

function [x,y,estr] = utm2ll(xi,yi,s_zone,i_type,ellipsoid) 

utmstruct = defaultm('utm'); 
utmstruct.zone = s_zone;        % e.g., '11S'

% if no ellipsoid is given then use Matlab default ellipsoid for that zone
if nargin == 4
    disp('no ellipsoid is provided as input -- using first listed one in utmgeoid.m');
    % NOTE: utmgeoid may return MULTIPLE geoids for a given zone
    [ellipsoid,estr] = utmgeoid(utmstruct.zone);
    %ellipsoid = almanac('earth','wgs84','meters');
    %ellipsoid = almanac('earth','clarke66','meters');
    %ellipsoid = [6.378206400000000e+06 0.082271854223002];

    ellipsoid = ellipsoid(1,:);
    estr0 = strtrim(estr(1,:));
    disp(sprintf('  ellipsoid %s : %.10e  %.10e',estr0,ellipsoid));
end
        
utmstruct.geoid = ellipsoid;        % assign ellipsoid
utmstruct = defaultm(utmstruct);
if i_type == 1
    [y,x] = minvtran(utmstruct,xi,yi);  % utm2ll
else
    [x,y] = mfwdtran(utmstruct,yi,xi);  % ll2utm
end

%------------------------------

if 0==1
    % to view UTM zones using Matlab GUI
    figure;
    axesm utm
    axesmui
    % THEN CLICK "ZONE"
    
    % to determine which zone a point is in, and the lat-lon bounds of the zone
    p1 = [32 -118];
    zs = utmzone(p1);
    [zonelats zonelons] = utmzone(zs);
    disp(sprintf('UTM zone %s extends from %.1f to %.1f in longitude and %.1f to %.1f in latitude',...
        zs,zonelons(1),zonelons(2),zonelats(1),zonelats(2)));
    
    % EXAMPLE
    format long, clc, clear, close all
    x0 = -118; y0 = 32;
    %x0 = 93.45; y0 = 7.65;
    p1 = [y0 x0];
    zs = utmzone(p1);
    [x1,y1] = utm2ll(x0,y0,zs,0);
    [x2,y2] = utm2ll(x1,y1,zs,1);
    disp(sprintf('original lon-lat point: %f, %f',x0,y0));
    disp(sprintf('recovered lon-lat point: %f, %f',x2,y2));
end

%==========================================================================
