%
% function ax1 = axes_expand(ax0,fac,iopt)
% Carl Tape, 16-Aug-2005
%
% This function inputs an axes box and outputs a new expanded axes box.
%
% EXAMPLE:
%    ax0 = [-121 -114]; ax1 = axes_expand(ax0,1.2,0)
%    ax0 = [-121 -114 31 37]; ax1 = axes_expand(ax0,1.2,0)
%    ax0 = [-121 -114 31 37]; ax1 = axes_expand(ax0,0.30,0)
%
% calls xxx
% called by xxx
%

function ax1 = axes_expand(ax0,fac,iopt)

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
    if iopt == 1, dinc = dy; else dinc = dx; end
    ax1(3) = ymin - dinc*(fac-1);
    ax1(4) = ymax + dinc*(fac-1);
    if ax1(4) <= ax1(3), ax1 = ax0; end
end
if ndim == 3
    zmin = ax0(3);
    zmax = ax0(4);
    dz = zmax-zmin;
    if iopt == 1, dinc = dz; else dinc = dx; end
    ax1(5) = zmin - dinc*(fac-1);
    ax1(6) = zmax + dinc*(fac-1);
    if ax1(6) <= ax1(5), ax1 = ax0; end
end

%===================================================
