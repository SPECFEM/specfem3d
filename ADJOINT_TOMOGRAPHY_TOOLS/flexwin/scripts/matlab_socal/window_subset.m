%
% window_subset.m
% CARL TAPE, 10-Jan-2008
% printed xxx
%
% This file reads in
%    (1) a list of event IDs and number of window picks
%    (2) another set of event IDs
% and it outputs
%    (1) a subset of all events and the window picks
%    (2) a subset of all events with at least Wmin window picks
%
% calls xxx
% called by xxx
%

clc
clear
close all

dir0 = '/net/sierra/raid1/carltape/results/WINDOWS/';

% read in list of window picks (window_summary.m)
%                          Z      R      T          NWIN
%    1    13935988 :     201    171    141   ==>     513
%    2    14155260 :     185    175    129   ==>     489
%file1 = [dir0 'T06_windows_sorted_by_event_pre_m0'];
%[index,eid0,junk1,Zwin,Rwin,Twin,junk2,winall] = textread(file1,'%f%s%s%f%f%f%s%f','headerlines',1);

file1 = [dir0 'T06_windows_sorted_by_event_pre_m0_matlab'];
[eid0,Zwin,Rwin,Twin,winall] = textread(file1,'%s%f%f%f%f');
win_mat = [Zwin Rwin Twin winall];

% read in list of event IDs
%file_eid = [dir0 'EIDs_fail'];
%file_eid = [dir0 'EIDs_pass'];
file_eid = '/net/sierra/raid1/carltape/results/SOURCES/socal_5/EIDs_all_loc_eid';

eids = textread(file_eid,'%s');
neid = length(eids);

%-----------------------------

% output file
ofile = [file_eid '_window'];
fid = fopen(ofile,'w');
for ii = 1:neid
    jj = strmatch(eids{ii},eid0,'exact');
    %disp([ii jj win_mat(jj,:)]);
    fprintf(fid,'%16s%6i%6i%6i%6i\n',char(eids(ii)),win_mat(jj,:));
end
fclose(fid);
disp([' writing ' ofile]);

% output file
Wmin = 500;      % KEY COMMAND
ofile = [file_eid '_window_' num2str(Wmin) '_plus'];
fid = fopen(ofile,'w');
for ii = 1:neid
    jj = strmatch(eids{ii},eid0,'exact');
    if win_mat(jj,end) >= Wmin
        fprintf(fid,'%16s%6i%6i%6i%6i\n',char(eids(ii)),win_mat(jj,:));
    end
end
fclose(fid);
disp([' writing ' ofile]);

%====================================================================
