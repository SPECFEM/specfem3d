%
% make_station_lists.m
% CARL TAPE, 20-Jan-2009
% printed xxx
%
% This program reads in a list of receivers and a list of events, and it
% outputs a set of lists sorted by azimuth and distance.  These are useful
% for sorting sets of data or plots, or also for a quick reference to find
% nearby stations for a particular event.
%
% calls xxx
% called by xxx
%

clc, clear, close all, format short, format compact

% output directory to dump the files into
odir = '/net/sierra/raid1/carltape/results/SOURCES/EID_STATION_LISTS/';

% format statement for output files
stfmt = '%12s%10.4f%10.4f%12s%10.4f%10.4f%10.4f%10.4f\n';

%------------------------------------------
% read in files
% Here are some example commands to make the input files:
%   awk '{print $2,$8,$9}' /net/sierra/raid1/carltape/results/SOURCES/socal_16/EIDs_lonlat_loc > EIDS_in
%   awk '{print $3"."$4,$1,$2}' /home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_gmt > STATIONS_ini

[eid_name,elon,elat] = textread('EIDS_in','%s%f%f');
[rec_name,rlon,rlat] = textread('STATIONS_in','%s%f%f');
neid = length(elon);
nrec = length(rlon);
disp(sprintf('%i sources and %i receivers',neid,nrec));

%------------------------------------------
% write files

for ii = 1:neid     % loop over events
    ename = eid_name{ii};
    elon0 = elon(ii);
    elat0 = elat(ii);
    disp(sprintf('Event %s : %i out of %i',ename,ii,neid));

    % distances and azimuths to all receivers
    [dists, azis] = distance(repmat(elat0,nrec,1),repmat(elon0,nrec,1),rlat,rlon);
    dists_km = dists*pi/180*6371;
    
    [junk1, idist] = sort(dists);
    [junk2, iaz] = sort(azis);

    % write files
    fid1 = fopen([odir 'STATIONS_by_dist_from_' ename],'w');
    fid2 = fopen([odir 'STATIONS_by_az_from_' ename],'w');
    for jj = 1:nrec
        k1 = idist(jj);
        k2 = iaz(jj);
        fprintf(fid1,stfmt,...
            rec_name{k1},rlon(k1),rlat(k1),ename,elon0,elat0,dists_km(k1),azis(k1));
        fprintf(fid2,stfmt,...
            rec_name{k2},rlon(k2),rlat(k2),ename,elon0,elat0,dists_km(k2),azis(k2));
    end
    fclose(fid1);
    fclose(fid2);
end

for jj = 1:nrec     % loop over receivers
    rname = rec_name{jj};
    rlon0 = rlon(jj);
    rlat0 = rlat(jj);
    disp(sprintf('Receiver %s : %i out of %i',rname,jj,nrec));

    % distances and (back-)azimuths to all events
    [dists, azis] = distance(repmat(rlat0,neid,1),repmat(rlon0,neid,1),elat,elon);
    dists_km = dists*pi/180*6371;
    
    [junk1, idist] = sort(dists);
    [junk2, iaz] = sort(azis);

    % write files
    fid1 = fopen([odir 'EIDS_by_dist_from_' rname],'w');
    fid2 = fopen([odir 'EIDS_by_az_from_' rname],'w');
    for ii = 1:neid
        k1 = idist(ii);
        k2 = iaz(ii);
        fprintf(fid1,stfmt,...
            eid_name{k1},elon(k1),elat(k1),rname,rlon0,rlat0,dists_km(k1),azis(k1));
        fprintf(fid2,stfmt,...
            eid_name{k2},elon(k2),elat(k2),rname,rlon0,rlat0,dists_km(k2),azis(k2));
    end
    fclose(fid1);
    fclose(fid2);
end

%====================================================================
