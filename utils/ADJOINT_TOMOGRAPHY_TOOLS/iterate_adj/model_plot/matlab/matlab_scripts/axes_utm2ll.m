%
% function ax_out = axes_utm2ll(ax_in,s_zone,i_type) 
% Carl Tape, 01-Aug-2007
%
% This function converts between a bounding region in lat-lon and in UTM.
%
% EXAMPLE:
%    ax_utm = axes_utm2ll([-121.6 -114.7 32.2 36.8],'11S',0) 
%
% calls utm2ll.m
% called by xxx
%

function ax_out = axes_utm2ll(ax_in,s_zone,i_type,ellipsoid) 

xmin0 = ax_in(1);
xmax0 = ax_in(2);
ymin0 = ax_in(3);
ymax0 = ax_in(4);

% if no ellipsoid is given then use Matlab default for that zone
if nargin == 3
    disp('no ellipsoid is provided as input');
    [ellipsoid,estr] = utmgeoid(s_zone)
end

[xmin,ymin] = utm2ll(xmin0,ymin0,s_zone,i_type,ellipsoid);
[xmax,ymax] = utm2ll(xmax0,ymax0,s_zone,i_type,ellipsoid);

ax_out = [xmin xmax ymin ymax];
        
if 0==1
    ax_ll = [-121.6 -114.7 32.2 36.8]
    ax_utm = axes_utm2ll(ax_ll,11,0) 
    axes_utm2ll(ax_utm,11,1) 
end

%==========================================================================
