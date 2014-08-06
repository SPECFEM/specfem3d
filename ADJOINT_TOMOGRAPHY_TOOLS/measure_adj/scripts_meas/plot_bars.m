%
% function plot_bars
% Carl Tape, 06-Nov-2007
%
% This plots a horizontal solid line, denoting a measurement,
   % and two dashed lines, denoting the uncertainty estimate.
% 
% calls xxx
% called by plot_mtm.m
%

function plot_bars(x1,x2,y1,sigma,cRGB,lsize)

xvec = [x1 x2];
yvec = [y1 y1];
plot(xvec,yvec,'linewidth',lsize,'color',cRGB)
plot(xvec,yvec+sigma,'--','linewidth',lsize/2,'color',cRGB)
plot(xvec,yvec-sigma,'--','linewidth',lsize/2,'color',cRGB)

%------------------------------
