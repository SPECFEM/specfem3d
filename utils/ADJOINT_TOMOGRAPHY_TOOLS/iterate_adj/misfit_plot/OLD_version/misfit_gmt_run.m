%
% misfit_gmt_run.m
% CARL TAPE, 03-Sept-2008
% printed xxx
%
%
% calls xxx
% called by xxx
%

%-------------------------

clear
close all
format compact

Tmin = 6;
Tmax = 40;
smodel1 = 'm00';
smodel2 = 'm07';
%stnm = 'USC'; netwk = 'CI';

ncomp = 3;
nper = 1;
nmodel = 2;

dir0 = pwd;
dir1 = [dir0 '/OUTPUT_MISFIT/EVENTS/'];

stations_file = [dir1 'STATIONS'];

%-------------------------
% read in list of event IDs and sources

nevent = 1;
eid = '14383980';
dir2 = [dir1 eid '/'];

%-----------------------
% read in stations file (used for synthetics)

[rlon0,rlat0,relev0,rburial0,stnm0,netwk0] = read_station_SPECFEM(stations_file);
nrec = length(stnm0);
%for ii = 1:nrec, strec{ii} = [stnm{ii} '.' netwk{ii}]; end
%strec = strec(:);

for irec = 1:nrec
   disp(sprintf('%i out of %i -- %s.%s',irec,nrec,stnm0{irec},netwk0{irec})); 
end

%-------------------------

nevent = 1; ievent = 1;
VR_array = zeros(nevent,nrec,ncomp,nper);

for irec = 1:nrec
    
    stnm = stnm0{irec};
    netwk = netwk0{irec};
    disp(sprintf('%i out of %i -- %s.%s',irec,nrec,stnm,netwk)); 
    
    odir = [dir2 stnm '.' netwk '/'];

    norms = misfit_gmt(odir,Tmin,Tmax,smodel1,smodel2,eid,stnm,netwk);

    VR_array(ievent,irec,1,1) = norms(end,1);  % Z
    VR_array(ievent,irec,2,1) = norms(end,2);  % R
    VR_array(ievent,irec,3,1) = norms(end,3);  % T
    
end

%=======================================================================
