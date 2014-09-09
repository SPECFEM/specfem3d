% 
% colors.m
% Carl Tape 01-Dec-2002
% printed xxx
%
% This gets the color-mapping matrices.
%
% calls cpt2cmap.m
% called by trigridN.m, plotcmapsN.m, raysAGU.m, xxx
%

% blue = (0,0,1)
% red = (1,0,0)
% green = (0,1,0)
% yellow = (0,1,1)
% white = (1,1,1)
% white = (0,0,0)

numc = 65;
white = ones(numc, 3);           % white (no good for contour plots)
bw = gray(numc);                 % black-to-white
wb = -gray(numc) + 1;            % white-to-black
nw = ceil(numc/2);
inc = 1/(nw-1);
br = [[ones(nw,1) [0:inc:1]' [0:inc:1]']; [[1-inc:-inc:0]' [1-inc:-inc:0]' ones(nw-1,1)]];
    % blue-to-red (top-to-bottom), with white in the middle
br = [[ones(nw,1) [0:inc:1]' [0:inc:1]']; [[1-inc:-inc:0]' [1-inc:-inc:0]' ones(nw-1,1)]];
    % blue-to-red (top-to-bottom), with white in the middle
for uy = 1:numc,
    rb(uy,:) = br(numc+1-uy, :);
end

% blue = (1,0,0)
% red = (0,0,1)
% green = (0,1,0)
% yellow = (0,1,1)
% white = (1,1,1)

% values from GMT 'seis' color palette
seis = [170	0	0;  206	0	0;  243	0	0;  255	24	0;  255	60	0;  255	97	0;
255	133	0; 255	170	0;  255	206	0;  255	243	0;  255	255	0; 255	255	0;
231	255	4;  161	255	17;  90	255	30;  51	249	64;  13	242	99;  0	194	152;
0	125	214;  0	68	248;  0	34	226];

% values for color palette used in TW phase velocity maps
ana = [202 80 0; 255 90 0; 255 110 0; 255 130 0; 255 150 0; 255 170 0; 255 190 0; ...
        255 205 0; 190 190 240; 170 170 240; 150 150 240; 130 130 240; 100 100 240; ...
        70 70 240; 30 30 220; 30 30 140];

br   = cpt2cmap(ana);
seis = cpt2cmap(seis);

%=========================================================================
