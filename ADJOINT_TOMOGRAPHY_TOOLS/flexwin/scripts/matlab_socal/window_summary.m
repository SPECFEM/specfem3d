%
% window_summary.m
% CARL TAPE, 03-Aug-2007
% printed xxx
%
% This file reads in a set of output files from the windowing code,
% compiles the statistics, and makes some plots.
%
% The user must provide a list of event IDs to read in.
%
% calls read_rectext.m, readCMT.m, readSTATIONS.m, socal_basemap.m
% called by xxx
%

clc
clear
close all

% the window tallies in the rectext file do NOT distinguish the channel
% type, but rather only the component (Z, R, T)
comps = {'Z','R','T'};
ncomp = length(comps);

%--------------------------------------
% USER INPUT

iwrite = 1;     % =1 to print out a list of events with too many windows
ifig = 0;       % =0 to avoid calling socal_basemap.m (which is not included here)

per = 6;        % CHANGE: shortest period used in the band-pass

stper = sprintf('T%2.2i',per);

% base directory
if 0==1
    dir0 = '/net/sierra/raid1/carltape/socal/socal_3D/TEST/';
    index = 1;
    dir = [dir0 stper '_' sprintf('%2.2i',index) '/'];
else
    dir = '/net/sierra/raid1/carltape/results/WINDOWS/model_m0/';
    %dir = '/net/sierra/raid1/carltape/results/WINDOWS/final_6s/';
    %dir = '/net/wagholi/scratch1/dchao/automeasure_work/collect_long_biot/socal/';
    %dir = '/net/wagholi/scratch1/dchao/automeasure_work/collect_long/socal/';
end

%--------------------------------------
% read in list of event IDs

dir_source = '/net/sierra/raid1/carltape/results/SOURCES/socal_5/';
eid_file = [dir_source 'SOCAL_FINAL_CMT_v5_eid'];
eid_all = textread(eid_file,'%s','headerlines',0);
nevent = length(eid_all);

% read in CMT solutions
cmt_file_all = [dir_source 'SOCAL_FINAL_CMT_v5'];
[date,tshift,hdur,slat,slon,dep,M,eid,elabel] = readCMT(cmt_file_all, 13, 1);

%--------------------------------------
% read in STATIONS file for synthetics
% This file shows ALL possible stations

stations_file = '/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem';
[rlon,rlat,relev,rburial,stnm,netwk] = readSTATIONS(stations_file);
nrec = length(stnm);
for ii = 1:nrec, strec{ii} = [stnm{ii} '.' netwk{ii}]; end
strec = strec(:);

% sort according to opt: 1 (station name), 2 (distance), 3 (azimuth)
% src_lon = slon(ii);
% src_lat = slat(ii);
% opt = 3;
% [rlon,rlat,relev,rburial,stnm,netwk,dist,az] = sortSTATIONS(stations_file,src_lon,src_lat,opt);
%writeSTATIONS(['qa' num2str(opt) ],rlon,rlat,relev,rburial,stnm,netwk,{0})

%--------------------------------------
% read in output files from windowing code

%nevent = 1;

% create indexing vector and array
ninds = nevent*nrec*ncomp;
index_array2vec = zeros(nevent,nrec,ncomp);
index_vec2array = zeros(ninds,3);

nw = 0;
for ievent = 1:nevent
    for irec = 1:nrec
        for icomp = 1:ncomp
            nw = nw+1;
            index_vec2array(nw,:) = [ievent irec icomp];
            index_array2vec(ievent,irec,icomp) = nw;
        end
    end 
end

% check that the operations are reversible
if 0==1
    index_vec2array( index_array2vec(3,23,2), : )
    a = index_vec2array(1214, :); index_array2vec(a(1),a(2),a(3))
    break
end

%---------------------

windows_array = zeros(nevent,nrec,ncomp);
ibmissing = zeros(nevent,1);

for ievent = 1:nevent       % loop over events 
    steid = eid_all{ievent};   
    src_lon = slon(ievent);
    src_lat = slat(ievent);
    disp([' Event ' num2str(ievent) ' out of ' num2str(nevent) ': ' steid]);
    
    if ifig==1
        figure;
        socal_basemap;
        plot(rlon,rlat,'.');
        text(rlon,rlat,stnm,'fontsize',7);
        plot(src_lon,src_lat,'p','markersize',20,'markeredgecolor','w','markerfacecolor','r');
        title([' Event ' steid '; depth ' sprintf('%.1f',dep(ievent)) ' km' ]);
    end
    
    ww = [steid '_T' sprintf('%2.2i',per) '_rectext_dist'];   % pre-sorted by arc-distance
    filename = [dir ww];
    
    if exist(filename,'file')
    
        % read in the number of windows picked
        [zrt_win,stnm1,netwk1] = read_rectext(filename);
        %win_tot_vec = sum(zrt_win,2);  % total window picks for each receiver
        %win_tot = sum(win_tot_vec);

        % list of receiver names
        nrec1 = length(stnm1);
        strec1 = [];
        for ii = 1:nrec1, strec1{ii} = [stnm1{ii} '.' netwk1{ii}]; end
        strec1 = strec1(:);

        % assign these window picks into windows_array
        for irec = 1:nrec
            % determine whether the receiver in the FULL FILE matches any
            % receiver in the windows file
            ind_rec = find( strcmp(strec(irec), strec1) == 1 );

            if ~isempty(ind_rec)
                for icomp = 1:ncomp
                    windows_array(ievent,irec,icomp) = zrt_win(ind_rec,icomp);
                end
            end
        end
        
    else
        ibmissing(ievent) = 1;
        disp([filename ' does not exist']);
    end
end

% indices of events with missing window picks
imissing = find(ibmissing==1);
ipresent = find(ibmissing==0);
nevent_present = length(ipresent);

%--------------------------------------------------------
% tabulations and figures

% sum by component
sum_comp = zeros(ncomp,1);
for icomp = 1:ncomp
    sum_comp(icomp) = sum(sum(windows_array(:,:,icomp)));
end

% sum by receiver
sum_rec = zeros(nrec,1);
for irec = 1:nrec
    sum_rec(irec) = sum(sum(windows_array(:,irec,:)));
end

% sum by event
sum_event = zeros(nevent,1);
for ievent = 1:nevent
    sum_event(ievent) = sum(sum(windows_array(ievent,:,:)));
end

% check that the sums are compatible
nwin = sum(sum_comp);
if length(unique([ sum(sum_comp) sum(sum_rec) sum(sum_event) ])) ~= 1
    error(' some problem with summing up windows_array');
else
    disp('  ');
    disp(sprintf(' Total number of windows (%i/%i events, %i receivers, %i components) : ',nevent_present,nevent,nrec,ncomp));
    disp(sprintf('%10s%10s%10s',comps{1},comps{2},comps{3}));
    disp(sprintf('%10i%10i%10i  ==> %10i',sum_comp,nwin));
end

% convert array to a data vector
windows_vec = zeros(ninds,1);
for ievent = 1:nevent
    for irec = 1:nrec
        for icomp = 1:ncomp
            ind = index_array2vec(ievent,irec,icomp);
            windows_vec(ind) = windows_array(ievent,irec,icomp);
        end
    end
end

% list records with greater than NWIN_MAX windows
nwin_max = 3;
imax = find(windows_vec >= nwin_max);
ofile = [stper '_records_with_' num2str(nwin_max) '_windows'];
disp('  '); disp([' writing ' ofile ' to file']);
fid = fopen(ofile,'w');
fprintf(fid,'%19s%8s%9s%6s%5s\n',' ','event','rec','comp','NWIN');
for ii=1:length(imax) 
    a = index_vec2array(imax(ii),:);
    ievent = a(1); irec = a(2); icomp = a(3);
    fprintf(fid,'%4i out of %4i : %8s%9s%6s%5i\n',...
        ii,length(imax),eid_all{ievent},strec{irec},comps{icomp},windows_vec(imax(ii)) );
end
fclose(fid);

% list window totals by event, sorted from most windows to least
sum_event_NaN = sum_event;
sum_event_NaN(imissing) = -1;
[junk, ievent_sort] = sort(sum_event_NaN,'descend');
ofile = [stper '_windows_sorted_by_event'];
disp([' writing ' ofile ' to file']);
fid = fopen(ofile,'w');
fprintf(fid,'%19s%7s%7s%7s%14s\n',' ',comps{1},comps{2},comps{3},'NWIN');
for ievent = 1:nevent
    ii = ievent_sort(ievent);
    fprintf(fid,'%4i%12s : %7i%7i%7i   ==> %7i\n',ievent,eid{ii},...
        sum(windows_array(ii,:,1)),sum(windows_array(ii,:,2)),...
        sum(windows_array(ii,:,3)),sum_event_NaN(ii) );
end
fprintf(fid,'%19s%7i%7i%7i%14i\n','TOTAL : ',sum_comp(1),sum_comp(2),sum_comp(3),nwin);
fclose(fid);

ofile = [stper '_windows_sorted_by_event_matlab'];
disp([' writing ' ofile ' to file']);
fid = fopen(ofile,'w');
for ievent = 1:nevent
    ii = ievent_sort(ievent);
    fprintf(fid,'%12s%7i%7i%7i%7i\n',eid{ii},...
        sum(windows_array(ii,:,1)),sum(windows_array(ii,:,2)),...
        sum(windows_array(ii,:,3)),sum_event_NaN(ii) );
end
fclose(fid);

% list window totals by receiver, sorted from most windows to least
[junk, irec_sort] = sort(sum_rec,'descend');
ofile = [stper '_windows_sorted_by_receiver'];
disp([' writing ' ofile ' to file']);
fid = fopen(ofile,'w');
fprintf(fid,'%18s%7s%7s%7s%14s\n',' ',comps{1},comps{2},comps{3},'NWIN');
for irec = 1:nrec
    ii = irec_sort(irec);
    fprintf(fid,'%5i%10s : %7i%7i%7i   ==> %7i\n',irec,strec{ii},...
        sum(windows_array(:,ii,1)),sum(windows_array(:,ii,2)),...
        sum(windows_array(:,ii,3)),sum_rec(ii) );
end
fprintf(fid,'%15s%7i%7i%7i%14i\n','TOTAL : ',sum_comp(1),sum_comp(2),sum_comp(3),nwin);
fclose(fid);

% list number of stations with >=1 window pick PER NETWORK
igood = find(sum_rec > 0);   % stations with at least ONE window pick
nets = unique(netwk(igood))
nunique = length(nets);
net_count_rec = zeros(nunique,1);
ofile = [stper '_network_totals'];
disp([' writing ' ofile ' to file']);
fid = fopen(ofile,'w');
for ii = 1:nunique
    net_count_rec(ii) = length(strmatch(nets{ii},netwk(igood),'exact'));
    fprintf(fid,'%4s : %3i\n',nets{ii},net_count_rec(ii));
end
fclose(fid);
 
% write a stations file showing all stations with at least ONE window pick
if 0==1
    filename = '/net/denali/home2/carltape/gmt/stations/seismic/Matlab_output/STATIONS_windows_m0';
    write_station_SPECFEM(filename,rlon(igood),rlat(igood),relev(igood),...
        rburial(igood),stnm(igood),netwk(igood));
end

% % plot histogram sorted by picks
% [sum_rec_sort, irec_sort] = sort(sum_rec,'descend');
% figure; hold on;
% plot(sum_rec_sort,'.');
% text([1:nrec],sum_rec_sort,stnm,'fontsize',6);
% grid on;
% xlabel(' Receiver'); ylabel([' Total number of windows (' num2str(nevent) ' events)']);

% display output for a particular event
if 0==1
    ievent = 133;
    for irec = 1:nrec
        disp([ sprintf('%6i : %8s %4i%4i%4i',irec,strec{irec},windows_array(ievent,irec,:)) ]);
    end
end

%====================================================================
