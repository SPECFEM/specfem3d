%
% readCMT.m
% CARL TAPE, 26-June-2007
% printed xxx
%
% This file reads a concatenated CMTSOLUTION file and outputs the data.
% The CMT tshift parameter is added to the origin time, which is the
% time used for the initial CMT solution.
%
% INPUT:
%   filename                 file containing all CMTSOLUTION files concatenated together
%   nlines_cmt               number of lines per CMT solution (13)
%   nspace_between_entries   number of spaces between CMTSOLUTION blocks (0 or 1)
%
% calls xxx
% called by slab_readCMT.m, socal_quakes_2007.m
%

function [date,tshift,hdur,lat,lon,dep,M,eid,elabel] = ...
    readCMT(filename, nlines_cmt, nspace_between_entries)

% read in concatenated CMTSOLUTION files
lines = textread(filename,'%s','delimiter','\n');
nlines = length(lines);

% number of events
enum = (nlines + nspace_between_entries) / (nlines_cmt + nspace_between_entries);
if mod(enum,1) ~= 0, error(' enum should be an integer'); end
disp([' File : ' filename ]);
disp([' Number of events : ' num2str(enum) ]);

% initialize vectors
date    = zeros(enum,1); tshift  = zeros(enum,1); hdur    = zeros(enum,1);
lat     = zeros(enum,1); lon     = zeros(enum,1); dep     = zeros(enum,1);
Mrr     = zeros(enum,1); Mtt     = zeros(enum,1); Mpp     = zeros(enum,1);
Mrt     = zeros(enum,1); Mrp     = zeros(enum,1); Mtp     = zeros(enum,1);

for kk=1:enum
    in1 = (kk-1)*(nlines_cmt + nspace_between_entries);
    
    % first full line of CMTSOLUTION, as a string
    ltemp = lines{in1+1};
    
    % replace W with space for more recent CMT solutions
    % SWEQ2006
    % PDEW2006
    if or(strcmp(ltemp(4),'W'),strcmp(ltemp(4),'Q')), ltemp(4) = ' '; end

    [j1,yr,mo,dy,hr,min,sec,lat_pde(kk),lon_pde(kk),dep_pde(kk),j5,j6,name1,name2,name3,name4,name5,name6] = ...
       strread(ltemp,'%s%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%s%s');
    [j1,j2,eid(kk)] = strread(lines{in1+2},'%s%s%s');
    [j1,j2,tshift(kk)] = strread(lines{in1+3},'%s%s%f');
    [j1,j2,hdur(kk)] = strread(lines{in1+4},'%s%s%f');
    [j1,lat(kk)] = strread(lines{in1+5},'%s%f');
    [j1,lon(kk)] = strread(lines{in1+6},'%s%f');
    [j1,dep(kk)] = strread(lines{in1+7},'%s%f');
    [j1,Mrr(kk)] = strread(lines{in1+8},'%s%f');
    [j1,Mtt(kk)] = strread(lines{in1+9},'%s%f');
    [j1,Mpp(kk)] = strread(lines{in1+10},'%s%f');
    [j1,Mrt(kk)] = strread(lines{in1+11},'%s%f');
    [j1,Mrp(kk)] = strread(lines{in1+12},'%s%f');
    [j1,Mtp(kk)] = strread(lines{in1+13},'%s%f');

    % NOTE the inclusion of the time shift here for the CMT solution
    date(kk) = datenum(yr,mo,dy,hr,min, sec + tshift(kk));
    
    elabel{kk} = [char(name1) ' ' char(name2) ' ' char(name3) ...
             ' ' char(name4) ' ' char(name5) ' ' char(name6)];
end

eid = eid(:);
elabel = elabel(:);

% moment tensor elements
M = [Mrr Mtt Mpp Mrt Mrp Mtp];

% compute moment based on Dahlen and Tromp, 1998, eq. 5.91
% compute moment magnitude based on Kanamori, 1977
% convention: r (up), theta (south), phi (east)
%M0 = 1/sqrt(2) * sqrt( Mrr.^2 + Mtt.^2 + Mpp.^2 + 2*(Mrt.*Mrt + Mrp.*Mrp + Mtp.*Mtp) );     
%Mw = (2/3)*log10(M0) - (2/3)*(11.8 - log10(5e-5));

%======================================================
