%
% function ax1 = axes_expand(ax0,fac)
% CARL TAPE, 16-Aug-2005
% printed xxx
%
% This function inputs an axes box and outputs a new axes box,
% expanded by the factor 'fac' in all 2,4,6 dimensions.
%
% EXAMPLE:
%    ax0 = [-121 -114]; ax1 = axes_expand(ax0,2)
%    ax0 = [-121 -114 31 37]; ax1 = axes_expand(ax0,2)
%    ax0 = [-121 -114 31 37]; ax1 = axes_expand(ax0,0.30)
%
% calls xxx
% called by xxx
%

function ax1 = axes_expand(ax0,fac)

% 1D, 2D, 3D
ndim = length(ax0)/2;
ax1 = zeros(1,ndim*2);

% return original axes if new axes are non-sensical

xmin = ax0(1);
xmax = ax0(2);
dx = xmax-xmin;
ax1(1) = xmin - dx*(fac-1);
ax1(2) = xmax + dx*(fac-1);
if ax1(2) <= ax1(1), ax1 = ax0; end

if ndim >= 2
    ymin = ax0(3);
    ymax = ax0(4);
    dy = ymax-ymin;
    ax1(3) = ymin - dy*(fac-1);
    ax1(4) = ymax + dy*(fac-1);
    if ax1(4) <= ax1(3), ax1 = ax0; end
end
if ndim == 3
    zmin = ax0(3);
    zmax = ax0(4);
    dz = zmax-zmin;
    ax1(5) = zmin - dz*(fac-1);
    ax1(6) = zmax + dz*(fac-1);
    if ax1(6) <= ax1(5), ax1 = ax0; end
end

%===================================================
