% 
% function x = seis_norm(t,y)
% Carl Tape, 03-Sept-2008
%
% This function computes the L2 norm of a seismogram.
% 
%             t2      2                  N        2      N        2
%          INT     s(t)   dt        dt  SUM  (s_i)      SUM  (s_i)
%   2         t1                        i=1             i=1
%  x   =    -------------------  = ----------------  = ---------
%               t2 - t1             dt*N                 N
%
% Note that the number will not depend on the length of the time series.
%
% calls xxx
% called by xxx
%

function x = seis_norm(t,y)

x = sqrt( sum(y.^2) / length(y) );

%============================================================
  