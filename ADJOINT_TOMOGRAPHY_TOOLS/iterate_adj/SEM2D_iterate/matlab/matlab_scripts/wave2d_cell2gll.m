%
% function 
% Carl Tape, 25-Jan-2010
%
% This function simply interpolates a function on a regular mesh (X,Y,F)
% onto a set of irregular points (xg,yg) to obtain fg.
%
% It is similar to wave2d_gll2cell.m, though NOT the opposite operation.
%
% calls xxx
% called by xxx
%

%function fg = wave2d_cell2gll(xg,yg,xc,yc,fc,nxc,nyc)
function fg = wave2d_cell2gll(xg,yg,X,Y,F)

% reshape xc, yc, fc into matrices
% NOTE: this assumes that they were ORIGINALLY made using meshgrid
%X = reshape(xc,nyc,nxc);
%Y = reshape(yc,nyc,nxc);
%F = reshape(fc,nyc,nxc);
fg = interp2(X,Y,F,xg,yg,'cubic');

%fg = interp2(X,Y,F,xg,yg,'cubic');

%=========================================================