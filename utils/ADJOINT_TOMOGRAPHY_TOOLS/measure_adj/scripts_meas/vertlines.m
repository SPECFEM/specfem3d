%
% function [xmat, ymat] = vertlines(xvec,ymin,ymax)
% CARL TAPE, 15-June-2005
% printed xxx
%
% This function inputs a vector of x-values and outputs
% two matrices that will plot vertical lines
% at the x-values specified in xvec.
%
% EXAMPLE:
%   [xmat,ymat] = vertlines(linspace(0,10,6),-1,4); figure; plot(xmat,ymat,'r');
%
% calls xxx
% called by LOTS OF FILES
% 

function [xmat, ymat] = vertlines(xvec,ymin,ymax)

n = length(xvec);
xvec = xvec(:)';        % xvec must be a row

% switch if incorrectLy entered
if ymin > ymax,
    temp = ymin;
    ymin = ymax;
    ymax = temp;
end

xmat = [xvec; xvec];
ymat = [ymin*ones(1,n); ymax*ones(1,n)];

%=====================================================
