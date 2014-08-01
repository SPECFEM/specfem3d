%
% function y = theoryHyp(m, param)
% CARL TAPE, 15-Nov-2003
% printed xxx
%
% Output is a set of y-data based on four
% parameters that describe the hyperbola.
%
% hyperbola: Ax + Bxy + Cy + D = 0
%
% Note that this is a LINEAR problem, Ax = b : 
%
%    |  x   x*y    y   1  | | A |   | 0 |
%    |  :    :     :   :  | | B | = | : |
%    |  :    :     :   :  | | C |   | : |
%    |  :    :     :   :  | | D |   | : |
%    |  :    :     :   :  |         | : |
%
% 
% See /matlab/scicomp/getellipse.m
%
% calls xxx
% called by c101D.m, testgenfitC.m, model_optimize2.m
%

function y = theoryHyp(m, param)

% input parameters (x-coordinates of data)
x = param(:);

% model parameters that we solve for
A = m(1); B = m(2); C = m(3); D = m(4);

% hyperbola: Ax + Bxy + Cy + D = 0
y = (-A*x - D) ./ (B*x + C);

%================================================
