%
% compare_misfit.m
% Carl Tape, 30-July-2009
%
% This program reads in two sets of window_chi files with IDENTICAL
% WINDOWS, computes the misfit for both, and then compute the variance
% reduction within each window.
%
% calls read_window_chi.m, plot_histo.m, readCMT.m, read_station_SPECFEM.m
% called by xxx
%

%-------------------------

clear
close all
format compact

% add path to additional matlab scripts
path(path,[pwd '/matlab_scripts']);

dir0 = '/home/carltape/RUNS/';

iVRlog = 1;    % variance reduction formula (1 for VRL; 0 for VR)
iwdiff = 1;    % demoninator for waveform difference (1 for new; 0 for old "standard" version)
wtag = sprintf('VR%i_wdiff%i',iVRlog,iwdiff);

iwrite = 0;
isciencehist = 0;    % 0 or 1
odirpub = '/home/carltape/manuscripts/2009/tomo_gji/latex/figures/scripts/hist_data/';

DTMIN = 1.0;   % see description below (1.0 or 0.0 for socal)

idataset = 3;        % 1 (tomo), 2 (extra), 3 (simulation)
dtags = {'tomo','extra','simulation'};
dtag = dtags{idataset};

%-------------------------

% min and max periods for the different bandpassed datasets
Tminvec = [2 3 6]; Tmaxvec = [30 30 30];
%Tminvec = [2]; Tmaxvec = [30];

comps = {'BHZ','BHR','BHT'};
    % regardless of the component label for the DATA, the measurement code
    % defaults so that the first two letters are always BH
ncomp = length(comps);
nper = length(Tminvec);
stBtag = 'BP';
for ii=1:nper
    stBtag = [stBtag num2str(Tminvec(ii))];
end

% strings for labeling
sTbp = repmat(cellstr(' '),1,nper);
for tt = 1:nper
    %sTbp{tt} = sprintf('[%is,%is]',Tminvec(tt),Tmaxvec(tt));
    sTbp{tt} = sprintf('%i-%is',Tminvec(tt),Tmaxvec(tt));
end
                
%-------------------------
% USER INPUT

%imod1 = input(' Enter the model 1 (0): ');
%imod2 = input(' Enter the model 2 (16): ');
imod1 = 0;
imod2 = 16;
stmod = sprintf('m%2.2i',imod2);     % put the VR files in the dir for the final model

idatacov = 2;
%idatacov = input(' Enter idatacov (1--by event, 2--by window, 3--none) : ');
    % 1: weight by N
    % 2: weight by E N_e in blocks that are N_e each in length
    % 3: weight by 1
    % -- N_e is the number of windows for event e
    % -- E is the number of events
    % -- N is the total number of windows: N = sum_e N_e
stdcovs = {'event','window','none'};
stdcov = stdcovs{idatacov};

ftag = [stmod '_' stdcov];
odir = ['OUTPUT_SUBSPACE/' stmod '/' stdcov '/'];

%iwrite = input(' Enter 1 to write files (0 otherwise) : ');

% files: CMTSOLUTIONS and stations
stsrcvers = '16';      % index for the set of sources (NOT the model iteration index)
dir_source = ['/home/carltape/results/SOURCES/socal_' stsrcvers '/'];
cmt_file_all = [dir_source 'SOCAL_FINAL_CMT_v' stsrcvers];
stations_file = '/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem';

% list of event IDs
%eid_file = [dir_source 'EIDs_only_loc'];
%eid_file = ['/net/sierra/raid1/carltape/results/EID_LISTS/syn_run_' stmod];
%eid_file = ['/net/sierra/raid1/carltape/results/EID_LISTS/kernels_run_' stmod];
%eid_file = ['/net/sierra/raid1/carltape/results/EID_LISTS/kernels_use_' stmod];

eid_file = ['/home/carltape/results/EID_LISTS/eids_' dtag];

%-------------------------
% read in list of event IDs and sources

% load the event IDs corresponding to the kernels
% load these as strings, since event IDs could have letters
eids = textread(eid_file,'%s');     
eids_num = textread(eid_file,'%f'); % numeric version
%eids_num = str2num(char(eids));
nevent = length(eids);

% read in CMT solutions
[date,tshift,hdur,slat,slon,dep,M,seid_cmt,elabel] = readCMT(cmt_file_all, 13, 1);
eid_cmt = str2num(char(seid_cmt));

% get the indices of these events into the full list
ievent_full = zeros(nevent,1);
for ii=1:nevent
    itemp = strmatch(eids(ii),seid_cmt);
    if ~isempty(itemp)
        ievent_full(ii) = itemp;
    else
        error(['event ' eids{ii} ' is not on the main list']);
    end
end
disp('all input EIDs matched to the main list');

%-----------------------
% read in stations file (used for synthetics)

[rlon,rlat,relev,rburial,stnm,netwk] = read_station_SPECFEM(stations_file);
nrec = length(stnm);
for ii = 1:nrec, strec{ii} = [stnm{ii} '.' netwk{ii}]; end
strec = strec(:);

% make list of networks
stnet = unique(netwk);
nnet = length(stnet);

%--------------------------------------
% create indexing array

% % create indexing vector and array
% ninds = nevent*nrec*ncomp*nper;
% index_array2vec = zeros(nevent,nrec,ncomp,nper);
% index_vec2array = zeros(ninds,4);
% 
% nw = 0;
% for ievent = 1:nevent
%     for irec = 1:nrec
%         for icomp = 1:ncomp
%             for iper = 1:nper
%                 nw = nw+1;
%                 index_vec2array(nw,:) = [ievent irec icomp iper];
%                 index_array2vec(ievent,irec,icomp,iper) = nw;
%             end
%         end
%     end 
% end
% nwin_tot = nw;    % = nevent*nrec*ncomp*nper

%-----------------------

iread = input(' Enter 1 to read in measurement files (0 otherwise) : ');
%iread = 1;
if iread == 1
    % range of events
    imin = 1; imax = nevent;        % default
    %imin = 1; imax = 10;
    %imin = 139; imax = imin;       
    %imin = 141; imax = imin;       % Chino Hills
    
    % model 1
    stmod = sprintf('m%2.2i',imod1);
    meas_array1 = read_window_chi_all(imin,imax,stmod,Tminvec,Tmaxvec,dir0,eids,strec,stnet,comps);
    save('meas_array1','meas_array1');

    % model 2
    stmod = sprintf('m%2.2i',imod2);
    meas_array2 = read_window_chi_all(imin,imax,stmod,Tminvec,Tmaxvec,dir0,eids,strec,stnet,comps);
    save('meas_array2','meas_array2');
else
   load(['meas_array1.mat']);
   load(['meas_array2.mat']);
end

% check at least that the station indices are identical
if sum( meas_array1(:,5) - meas_array2(:,5) ) ~= 0
    error('mismatch of the two window files');
end

% total number of windows
N = length(meas_array1);

% combine the columns that you want into one datafile
% GET RID OF THE 0.5 FACTOR IN tr_chi, SO THAT IT IS DT^2 / SIGMA^2
meas_all = [meas_array1(:,[1:7 9 10 13 14 11]) ...     % first 11 columns
    meas_array1(:,12) meas_array2(:,12) meas_array1(:,8) meas_array2(:,8) ...
    meas_array1(:,15) meas_array2(:,15) meas_array1(:,16) meas_array2(:,16) ...
    meas_array1(:,17) meas_array2(:,17) meas_array1(:,18) meas_array2(:,18) ...
    meas_array1(:,19) meas_array2(:,19) meas_array1(:,23) meas_array2(:,23) ...
    meas_array1(:,25) meas_array2(:,25)];
%  1  kinds
%  2  index_T
%  3  index_event
%  4  index_network
%  5  index_rec
%  6  index_comp
%  7  iwin
%  8  seisdur
%  9  windur
% 10  seisd2
% 11  wind2
% 12  T_pmax_dat
% 13  m00 - T_pmax_syn
% 14  m16 - T_pmax_syn
% 15  m00 - iker
% 16  m16 - iker
%-------------------
% waveform difference information
% 17  m00 - seiss2
% 18  m16 - seiss2
% 19  m00 - wins2
% 20  m16 - wins2
% 21, m00 - seis_diff2
% 22, m16 - seis_diff2
% 23, m00 - win_diff2
% 24, m16 - win_diff2
%-------------------
% 25, m00 - DT-CC
% 26, m16 - DT-CC
% 27, m00 - DT-MT
% 28, m16 - DT-MT
% 29, m00 - tr_chi
% 30, m16 - tr_chi
%-------------------
% 31, m00 - seis_chi
% 32, m16 - seis_chi
% 33, m00 - win_chi
% 34, m16 - win_chi
% 35 --> VRseis wdiff (with NaN)
% 36 --> VRwin wdiff (with NaN)
%-------------------

% find the unique seismograms (value from first window in each record)
bseis1 = (meas_all(:,7) == 1);
nseis = sum(bseis1);
disp(sprintf('%i out of %i windows are the only window on the seismogram',nseis,N));

for kk = 0:1    % loop over both models
    % find the number of MT and CC measurements used, and replace MT-DT with the CC value
    iMT = find(meas_all(:,15+kk) == 1);
    iCC = find(meas_all(:,15+kk) == 2);
    nMT = length(iMT);
    disp(sprintf('%i out of %i windows use multitaper measurment',nMT,N));

    % compare DT for MT vs CC
    figure; hold on; ymx = 4;
    plot(meas_all(:,25+kk),meas_all(:,27+kk),'.');
    plot(meas_all(iCC,25+kk),meas_all(iCC,27+kk),'ro');
    plot([-1 1]*ymx,[-1 1]*ymx,'r--'); grid on; axis equal, axis([-1 1 -1 1]*ymx);
    xlabel('Cross-correlation DT'); ylabel('Multitaper DT (when possible)');
    
    % KEY: replace MT DT zero values with the CC DT values
    meas_all(iCC,27+kk) = meas_all(iCC,25+kk);
    plot(meas_all(iCC,25+kk),meas_all(iCC,27+kk),'ko');
end

%itemp = find(and(meas_all(:,23)==-2.9,meas_all(:,25) > 0))
%display_meas(meas_array2(itemp,:),Tminvec,eids,strec,comps);
%break

%----------------------------------------
% WAVEFORM DIFFERENCE MISFIT FUNCTION

if iwdiff == 1
    % new version
    seis_chi1 = meas_all(:,21) ./ sqrt( meas_all(:,10) .* meas_all(:,17) );
    seis_chi2 = meas_all(:,22) ./ sqrt( meas_all(:,10) .* meas_all(:,18) );
    win_chi1  = meas_all(:,23) ./ sqrt( meas_all(:,11) .* meas_all(:,17) );
    win_chi2  = meas_all(:,24) ./ sqrt( meas_all(:,11) .* meas_all(:,18) );  
else
    % old "standard" version
    seis_chi1 = meas_all(:,21) ./ meas_all(:,10);
    seis_chi2 = meas_all(:,22) ./ meas_all(:,10);
    win_chi1  = meas_all(:,23) ./ meas_all(:,11);
    win_chi2  = meas_all(:,24) ./ meas_all(:,11);  
end

% add misfit to meas_all
meas_all = [meas_all seis_chi1 seis_chi2 win_chi1 win_chi2];

%----------------------------------------
% VARIANCE REDUCTION

if iVRlog == 1
    VRseis = log( seis_chi1 ./ seis_chi2 );
    VRwin = log( win_chi1 ./ win_chi2 );
    edges_vr = [-4:0.5:4]; ylims = [0 0.4];
    VRtag = 'VRL';
else
    VRseis = (seis_chi1 - seis_chi2) ./ seis_chi1;
    VRwin = (win_chi1 - win_chi2) ./ win_chi1;   
    edges_vr = [-4:0.25:1]; ylims = [0 0.4];
    VRtag = 'VR';
end
    

% Remove windows that are have time shifts less than DTMIN both in the
% initial model and in the final model -- we are not interested in the VR
% for these records.  But we are interested if the time shift gets WORSE
% than DTMIN in going from m00 to m16.
%DTMIN = 0.0;
DTMIN_tag = sprintf('DTMIN_%3.3ims',DTMIN*1000);
bbad = and( abs(meas_all(:,25)) < DTMIN, abs(meas_all(:,26)) < DTMIN );
bgood = ~bbad;
ngood = sum(bgood);
VRwin(bbad) = NaN;

if DTMIN == 0
    VRseis(~bseis1) = NaN;    % set VR for non-unique entries to NaN
else
    disp('excluding seismograms that have windows fit before and after');
    if 0==1
        % exclude seismograms in which ANY window has a time shifts less than DTMIN
        % both in the initial model and in the final model
        VRseis(~bseis1) = NaN;    % set VR for non-unique entries to NaN
        VRseis(bbad) = NaN;      % set VR to NaN if ANY window is NaN
    else
        % exclude seismograms in which ALL windows have time shifts less than DTMIN
        % both in the initial model and in the final model -- THIS IS SLOW

        %for ii = 1:N        % loop over all windows
        ind1=3; ind2=5; ind3=6; ind4=2;
        for ii = 1:N
        %for ii = 10:10
           if mod(ii,1000) == 0, disp(sprintf('%i out of %i',ii,N)); end
           sind = meas_all(ii,[ind1 ind2 ind3 ind4]);    % index into record

           % find the other windows on the seismogram and assign NaN if all windows
           % on the seismogram have VRwin NaN
           imatch = find( and( and(meas_all(:,ind1)==sind(1),meas_all(:,ind2)==sind(2)),...
               and(meas_all(:,ind3)==sind(3),meas_all(:,ind4)==sind(4)) ) == 1);
           inans = isnan(VRwin(imatch));
           if all(inans)                    % all windows are bad
               VRseis(imatch) = NaN;
           else                             % at least one window is good
               ikeeper = imatch(~inans);
               VRseis(setdiff(imatch,ikeeper(1))) = NaN;
           end
        end
    end
end

%temp = [VRseis VRwin]; disp(temp(1:30,:)), sum(~isnan(VRseis)), sum(~isnan(VRwin))

% good seismograms, unrepeated
bgoodseis = ~isnan(VRseis);
ngoodseis = sum(bgoodseis);

disp(sprintf('%i out of %i windows have DT > %.2f s for m%2.2i or m%2.2i',ngood,N,DTMIN,imod1,imod2));
disp(sprintf('%i out of %i seismograms have DT > %.2f s for m%2.2i or m%2.2i',ngoodseis,nseis,DTMIN,imod1,imod2));

% Remove windows that are worse within records that are better
%bbad2 = and( VRseis > 0, VRwin <  0); length(bbad2)
%display_meas(meas_array1(bbad2==1,:),Tminvec,eids,strec,comps);
%break

% fraction of records and windows that have a worse waveform difference
fworse_seis = length(find(VRseis(bgoodseis) < 0)) / ngoodseis;
fworse_win = length(find(VRwin(bgood) < 0)) / ngood;

% add VR to meas_all
meas_all = [meas_all VRseis VRwin];

VRvals = [mean(VRwin(bgood)) median(VRwin(bgood)) mean(VRseis(bgoodseis)) median(VRseis(bgoodseis))];
disp(sprintf('Variance reduction values for the %i windows:',ngood));
disp(sprintf('%12s%12s%12s%12s','VRwin-mean','VRwin-med','VRseis-mean','VRseis-med'));
disp(sprintf('%12.4f%12.4f%12.4f%12.4f',VRvals));

figure; nr=3; nc=2;
edges_chi = [0:0.25:4];

subplot(nr,nc,1);
x = win_chi1(bgood);
%x = log(win_chi1(bgood));
Nbin = plot_histo(x,edges_chi);
xlabel('chi2 for each window'); ylim(ylims); grid on;
title(sprintf('chi2 for m%2.2i WINDOWS: mean %.3f, median %.3f',imod1,mean(x),median(x)));

subplot(nr,nc,2);
x = seis_chi1(bgoodseis);
%x = log(seis_chi1(bgoodseis));
Nbin = plot_histo(x,edges_chi);
xlabel('chi2 for each seismogram'); ylim(ylims); grid on;
title(sprintf('chi2 for m%2.2i SEIS: mean %.3f median %.3f',imod1,mean(x),median(x)));

subplot(nr,nc,3);
x = win_chi2(bgood);
%x = log(win_chi2(bgood));
Nbin = plot_histo(x,edges_chi);
xlabel('chi2 for each window'); ylim(ylims); grid on;
title(sprintf('chi2 for m%2.2i WIN: mean %.3f median %.3f',imod2,mean(x),median(x)));

subplot(nr,nc,4);
x = seis_chi2(bgoodseis);
%x = log(seis_chi2(bgoodseis));
Nbin = plot_histo(x,edges_chi);
xlabel('chi2 for each seismogram'); ylim(ylims); grid on;
title(sprintf('chi2 for m%2.2i SEIS: mean %.3f median %.3f',imod2,mean(x),median(x)));

subplot(nr,nc,5);
Nbin = plot_histo(VRwin(bgood),edges_vr);
xlabel([VRtag ' for each window']); ylim(ylims); grid on;
title(sprintf('WIN: %.2f > 0, %smean = %.3f, %smed = %.3f',1-fworse_win,VRtag,VRvals(1),VRtag,VRvals(2)));

subplot(nr,nc,6);
Nbin = plot_histo(VRseis(bgoodseis),edges_vr);
xlabel([VRtag ' for each seismogram']); ylim(ylims); grid on;
title(sprintf('SEIS: %.2f > 0, %smean = %.3f, %smed = %.3f',1-fworse_seis,VRtag,VRvals(3),VRtag,VRvals(4)));

fontsize(10), orient tall, wysiwyg

%---------------------
% MT-DT and CC-DT for final model

DTvals = zeros(1,4);
x1 = meas_all(bgoodseis,26);
x2 = meas_all(bgoodseis,28);
DTvals = [mean(x1) std(x1) mean(x2) std(x2)];

if 1==1
    figure; nr=2; nc=1; ylims = [0 0.3];
    edges_DT = [-3:0.25:3];
    tlabs = {'CC','MT'};
    
    for kk=[0 2]
        % windows
        subplot(nr,nc,1+kk/2);
        x = meas_all(bgood,26+kk);
        Nbin = plot_histo(x,edges_DT);
        xlabel('DT for each window'); ylim(ylims); grid on;
        title(sprintf('%s-DT, m%2.2i, %i events: mean %.3f pm %.3f',tlabs{1+kk/2},imod2,nevent,mean(x),std(x)));

%         % seismograms
%         subplot(nr,nc,2+kk);
%         x = meas_all(bgoodseis,26+kk);
%         Nbin = plot_histo(x,edges_DT);
%         xlabel('DT for each window'); ylim(ylims); grid on;
%         title(sprintf('%s-DT, m%2.2i, %i events: mean %.3f pm %.3f',tlabs{1+kk/2},imod2,nevent,mean(x),std(x)));
    end
    fontsize(10), orient tall, wysiwyg
end

% %----------------------------------------------
% % check some of the output window rows
% if 0==1
%     itestvec = [41];
%     for ii = 1:length(itestvec)
%         disp('====>');
%         itest = itestvec(ii);
%         meas_array1(itest,14)
%         meas_array2(itest,14)
%         display_meas(meas_array1(itest,:),Tminvec,eids,strec,comps);
%         display_meas(meas_array2(itest,:),Tminvec,eids,strec,comps);
%     end
% end
  
% add a couple more measures
% 32 - Pwin, power in window relative to record
% 33 - Pwin*VRwin, quality of VR for window (will favor single-window records)
%meas_all_mod = [meas_all meas_all(:,11)./meas_all(:,10) meas_all(:,11)./meas_all(:,10).*meas_all(:,31)];
meas_all_mod = meas_all;

% if isciencehist==1
% 
%     % 1  - ii
%     % 2  - eid
%     % 3  - ind
%     % 4  - sta
%     % 5  - net
%     % 6  - comp
%     % 7  - iwin
%     % 8  - m16 - DT-CC
%     % 9  - m16 - DT-MT
%     % 10 - m00 - seis_chi
%     % 11 - m16 - seis_chi
%     % 12 - m00 - win_chi
%     % 13 - m16 - win_chi
% 
%     %dsort = meas_all_mod(bgood,:);        % extract a subset
%     dsort = meas_all_mod;
%     klabs = {'win','seis'};
%     for kk = 1:2
%         fpre = 'hist_all_wdiff';
%         %fpre = 'hist_T006_T030_DT';
%         
%         fname = [odirpub fpre '_'  DTMIN_tag '_' dtag '_' wtag '_' klabs{kk} '.dat'];
%         fid = fopen(fname,'w');
%         pp = 0;
%         for ii = 1:length(dsort)
%             ilist = 0;
%             if any(kk==[2 4])          % seismogram list
%                 %if and(~isnan(dsort(ii,[33 34])),dsort(ii,7)==1), ilist = 1; end
%                 if ~isnan(dsort(ii,35)), ilist = 1; end
%             else                       % window list
%                 if ~isnan(dsort(ii,36)), ilist = 1; end
%             end
%                 
%             if ilist==1
%                 pp = pp+1;
%                 fprintf(fid,['%6i%12i%6i%6s%6s%6s%6i' repmat('%8.3f',1,6) '\n'],...
%                     pp,eids_num(dsort(ii,3)),dsort(ii,1),...
%                     char(stnm{dsort(ii,5)}),char(stnet{dsort(ii,4)}),char(comps{dsort(ii,6)}),...
%                     dsort(ii,[7 26 28 29:32]));
%             end
%         end
%         fclose(fid);
%     end
% end

if iwrite == 1
    
    if isciencehist==0
        klabs = {'DT2','DT_chi2','seis_chi1','seis_chi2','win_chi1','win_chi2','VRseis','VRwin'};
        kind = [-26 30 31 32 33 34 -35 -36];          % columns to sort, corresponding to klabs
    else
        odir = odirpub;
        klabs = {'seis_chi1','win_chi1'};
        kind = [31 33];
    end
    numk = length(klabs);
    
    disp('writing out sorted files...');
    for kk = 1:numk
    %for kk = numk:numk
    %for kk = 1:1
        dsort = sortrows(meas_all_mod,kind(kk));
        fid = fopen([odir ftag '_' stBtag '_' DTMIN_tag '_' wtag '_' dtag '_sort_' klabs{kk} '.txt'],'w');
        fprintf(fid,['%6s%12s%8s%6s%6s%6s%10s%6s' repmat('%8s',1,12) '\n'],' ','eid','ind',...
            'stnm','net','comp','BP','iwin','DTCC1','DTCC2','DTMT1','DTMT2',...
            'DTchi1','DTchi2','Fseis1','Fseis2','Fwin1','Fwin2','VRseis','VRwin');
        pp = 0;
        for ii = 1:N
            ilist = 0;
            if any(abs(kind(kk))==[31 32 35])          % seismogram list
                %if and(~isnan(dsort(ii,[33 34])),dsort(ii,7)==1), ilist = 1; end
                if ~isnan(dsort(ii,35)), ilist = 1; end
            else                       % window list
                if ~isnan(dsort(ii,36)), ilist = 1; end
            end
            if ilist == 1
                pp = pp+1;
                fprintf(fid,['%6i%12i%8i%6s%6s%6s%10s%6i' repmat('%8.3f',1,12) '\n'],...
                    pp,eids_num(dsort(ii,3)),...
                    dsort(ii,1),...
                    char(stnm{dsort(ii,5)}),...
                    char(stnet{dsort(ii,4)}),...
                    char(comps{dsort(ii,6)}),...
                    char(sTbp{dsort(ii,2)}),...
                    dsort(ii,[7 25:36]));
            end
        end
        %fprintf(fid,'%12s%8i%6s%6s%6s%12.4f%12s%12s\n','TOTAL',sum(Ns),'--','--','--',sum(dnorm_sq),'--','--');
        fclose(fid);
    end
    
    % values to go with the histograms
    % 1 - m16 DTCC mean
    % 2 - m16 DTCC std
    % 3 - m16 DTMT mean
    % 4 - m16 DTMT std
    % 5 - VRwin mean
    % 6 - VRwin median
    % 7 - VRseis mean
    % 8 - VRseis median
    vtemp = [DTvals VRvals];
    disp('writing out statistic values...');
    fid = fopen([odir ftag '_' stBtag '_' DTMIN_tag '_' wtag '_' dtag '_values.txt'],'w');
    fprintf(fid,[repmat('%12.3e',1,length(vtemp)) '\n'],vtemp);
    fclose(fid);
end

%=======================================================================
