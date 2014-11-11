%
% function outmat = cpt2cmap(inmat)
% Carl Tape, 14-March-2005
% printed xxx
%
% This function inputs a 3 column matrix of colors (RGB)
% and outputs a matrix suitable for Matlab plots.
%
% calls xxx
% called by xxx
%

function outmat = cpt2cmap(inmat)

numc = 65;
len = length(inmat);
temp = inmat/max(max(inmat));  % normalize to 0-1 range
rv = temp(:,1);
gv = temp(:,2);
bv = temp(:,3);

% I don't know a better way to do this -- this interpolates to numc points.
numc = 65;
outmat = [  interp1(linspace(1,numc,len), rv', [1:numc])' ...
            interp1(linspace(1,numc,len), gv', [1:numc])' ...
            interp1(linspace(1,numc,len), bv', [1:numc])'];

%==================================================
