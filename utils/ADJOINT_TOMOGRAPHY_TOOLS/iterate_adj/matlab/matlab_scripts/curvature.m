% 
% function [i0,kap] = curvature(x,y)
% CARL TAPE, 25-July-2005
%
% This function compute the curvature of a 1-D function
% and returns the computed curvature, as well as the
% point of maximum POSITIVE curvature.
%
% This was written for quantifying L-curves.
%
% calls xxx
% called by spline_wang_D.m, test_del2.m
%

function [i0,kap] = curvature(x,y)

f1 = gradient(y,x);
f2 = gradient(f1,x);
kap = f2 ./ (1 + f1.^2).^(3/2);  % see Mathworld

%[kap0,i0] = max(abs(kap));
[kap0,i0] = max(kap);

%============================================================
  