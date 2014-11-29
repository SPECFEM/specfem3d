%
% function [xbox,ybox] = boxpts(ax0,n)
% CARL TAPE, 06-Sept-2005
% printed xxx
%
% This function inputs an axes box and outputs a set of points describing
% the boundary of the box.
%
% EXAMPLE:
%   ax0 = [-121.6 -114.7 32.2 36.8]; [xbox,ybox] = boxpts(ax0,100); figure, plot(xbox,ybox,'.');
%
% calls xxx
% called by xxx
%

function [xbox,ybox] = boxpts(ax0,n)

xmin = ax0(1);
xmax = ax0(2);
ymin = ax0(3);
ymax = ax0(4);

% total length of box
xran = xmax - xmin;
yran = ymax - ymin;
len = 2*xran + 2*yran;

% spacing increment
dx = len/n;
nx = round(xran/dx);
ny = round(yran/dx);
xpts = linspace(xmin,xmax,nx);
ypts = linspace(ymin,ymax,ny);

xbox = [xpts ones(1,ny)*xmax fliplr(xpts) ones(1,ny)*xmin ]';
ybox = [ones(1,nx)*ymin ypts ones(1,nx)*ymax fliplr(ypts) ]';

%===================================================
