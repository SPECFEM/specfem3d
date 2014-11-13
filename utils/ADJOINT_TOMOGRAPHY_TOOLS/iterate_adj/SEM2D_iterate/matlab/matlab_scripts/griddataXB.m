%
% function [X,Y,Z] = griddataXB(xvec,yvec,zvec,npts,stype)
% Carl Tape, 09-July-2004
% printed xxx
%
% Copied from griddataX.m on 08-May-2004
% Converts irregular points to regular points, for mesh plotting.
%
% INPUT:
%   xvec,yvec   = map coordinates of irregular points
%   zvec        = function values at (xi,yi)
%   npts        = number of points in horizontal interpolation
%   stype       = type of interpolation
%
% 		'linear'    - Triangle-based linear interpolation (default).
% 		'cubic'     - Triangle-based cubic interpolation.
% 		'nearest'   - Nearest neighbor interpolation.
% 		'v4'        - MATLAB 4 griddata method.
%
% OUTPUT:
%   X,Y         = interpolated mesh
%   Z           = interpolated function
%
% calls xxx
% called by c164E.m, GPS_A.m, omsplotsF.m
%

function [X,Y,Z] = griddataXB(xvec, yvec, zvec, npts, stype)

% construct mesh with UNIFORM spacing in x and y directions
xlin  = linspace(min(xvec), max(xvec), npts);
dx    = xlin(2) - xlin(1);
ylin  = [min(yvec) : dx : max(yvec)+dx];
[X,Y] = meshgrid(xlin,ylin);

%zvec = zvec(:)';    % row vector

% determine interpolated function using xvec,yvec input
Z = griddata(xvec, yvec, zvec, X, Y, stype);

%==================================================
