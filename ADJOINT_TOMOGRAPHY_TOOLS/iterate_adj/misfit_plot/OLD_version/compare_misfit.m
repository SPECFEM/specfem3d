%
% compare_misfit.m
% CARL TAPE, 23-July-2008
% printed xxx
%
% This script reads in a set of window_chi output measurement files from
% mt_measure_adj.f90, and it tabulates and plots misfit function values.
%
% /cig/seismo/3D/ADJOINT_TOMO/iterate_adj/matlab/
%
% Copied from compute_misfit.m on 23-July-2008
%
% calls read_window_chi.m, plot_histo.m, readCMT.m, read_station_SPECFEM.m
% called by xxx
%

%-------------------------

clear
close all
format compact

dir0 = '/net/sierra/raid1/carltape/socal/socal_3D/RUNS/';

Tvec = [6 2];
    % LOWEST period used for each measurement band-pass (upper is 40 s)
comps = {'BHZ','BHR','BHT'};
    % regardless of the component label for the DATA, the measurement code
    % defaults so that the first two letters are always BH
ncomp = length(comps);
nper = length(Tvec);

%DT = 0.011;     % really this is only needed because it was left out from
                % integrating the waveforms in mt-measure_adj.f90

%-------------------------
% USER INPUT

iwrite = input(' Enter 1 to write files (0 otherwise) : ');

% files: event IDs, CMTSOLUTIONS, stations
dir_source = '/net/sierra/raid1/carltape/results/SOURCES/socal_09/';
%eid_file = [dir_source 'EIDs_only_loc'];
%eid_file = ['/net/sierra/raid1/carltape/results/EID_LISTS/kernels_use_' stmod];
eid_file = '/net/sierra/raid1/carltape/results/EID_LISTS/kernels_use_m05';
cmt_file_all = [dir_source 'SOCAL_FINAL_CMT_v09'];
stations_file = '/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem';

%-------------------------
% read in list of event IDs and sources

% load the event IDs corresponding to the kernels
% load these as strings, since event IDs could have letters
eids = textread(eid_file,'%s');     
eids_num = textread(eid_file,'%f'); % numeric version
%eids_num = str2num(char(eids));
nevent = length(eids);

% read in CMT solutions
[date,tshift,hdur,slat,slon,dep,M,eid,elabel] = readCMT(cmt_file_all, 13, 1);

%-----------------------
% read in stations file (used for synthetics)

[rlon,rlat,relev,rburial,stnm,netwk] = read_station_SPECFEM(stations_file);
nrec = length(stnm);
for ii = 1:nrec, strec{ii} = [stnm{ii} '.' netwk{ii}]; end
strec = strec(:);

% make list of networks
stnet = unique(netwk);
nnet = length(stnet);

%-----------------------

iread = input(' Enter 1 to read in measurement files (0 otherwise) : ');

imodel_min = input(' Enter minimum model index (0-5) : ');
imodel_max = input(' Enter maximum model index (1-5) : ');
nmodel = imodel_max - imodel_min;

wdiff_array = zeros(nevent,nrec,ncomp,nper,nmodel);

if iread == 1

    % range of events
    imin = 1; imax = nevent;        % default
    %imin = 4; imax = 7;
    %imin = 4; imax = imin;
    
    for imod = imodel_min:imodel_max
        
        stmod = sprintf('m%2.2i',imod);
        disp('====================================');
        disp(['MODEL ' stmod]);
    
        k2 = 0;
        meas_array = [];
        %for tt = 1:1
        for tt = 1:nper
            Tper = Tvec(tt);
            sTper = sprintf('%3.3i',Tper);
            if imod <= 2, sTper2 = sprintf('%2.2i',Tper); else sTper2 = sTper; end   % CHT
            disp('-----------------------------------------');

            for ii = imin:imax    % loop over events

                disp(sprintf('Event %i out of %i -- %s',ii,nevent,eids{ii}));
                dir1 = [dir0 eids{ii} '/' stmod '/MEASURE_T' sTper '/'];
                mfile = [eids{ii} '_T' sTper2 '_window_chi'];
                filename = [dir1 mfile];

                if ~exist(filename,'file')
                    disp('--> file does not exist (or nwin = 0)');
                else
                    % read the window_chi file
                    [stnet0,strec0,comp,iwin,iker,t1,t2,...
                    chiMT_dT,chiMT_dA,chiCC_dT,chiCC_dA,...
                    measMT_dT,measMT_dA,measCC_dT,measCC_dA,...
                    sigmaMT_dT,sigmaMT_dA,sigmaCC_dT,sigmaCC_dA,...
                    wind2,wins2,win_chi0,windur, recd2,recs2,rec_chi0,recdur,...
                    tr_chi,am_chi] = read_window_chi(filename);

                    % waveform differences normalized by data-squared
                    win_chi = win_chi0 ./ wind2;
                    rec_chi = rec_chi0 ./ recd2;

                    nwin = length(strec0);
                    disp(['--> nwin = ' num2str(nwin)]);

                    % assign the string network an integer index
                    index_net = NaN * ones(nwin,1);
                    for inet = 1:nnet
                        ind_net = find( strcmp(stnet(inet), stnet0) == 1 );
                        index_net(ind_net) = inet;
                    end

                    % assign the string receiver an integer index
                    index_rec = NaN * ones(nwin,1);
                    for irec = 1:nrec
                        ind_rec = find( strcmp(strec(irec), strec0) == 1 );
                        index_rec(ind_rec) = irec;
                    end

                    % assign the string component an integer index
                    index_comp = NaN * ones(nwin,1);
                    for icomp = 1:ncomp
                        ind_comp = find( strcmp(comps(icomp), comp) == 1 );
                        index_comp(ind_comp) = icomp;
                    end

                    % measurement indices
                    k1 = k2+1;
                    k2 = k1+nwin-1;
                    kinds = [k1:k2]';

                    meas_array(kinds,:) = [kinds tt*ones(nwin,1) ii*ones(nwin,1) index_net index_rec index_comp iwin ...
                        iker measCC_dT sigmaCC_dT measCC_dA sigmaCC_dA  measMT_dT sigmaMT_dT win_chi rec_chi tr_chi];
                    %  1  kinds
                    %  2  index_T
                    %  3  index_event
                    %  4  index_network
                    %  5  index_rec
                    %  6  index_comp
                    %  7  iwin
                    %  8  iker
                    %  9  measCC_dT
                    % 10  sigmaCC_dT
                    % 11  measCC_dA
                    % 12  sigmaCC_dA
                    % 13  measMT_dT
                    % 14  sigmaMT_dT
                    % 15  win_chi
                    % 16  rec_chi                
                    % 17  tr_chi
                end

            end
        end

        for tt = 1:nper
            for ii = imin:imax
                disp(sprintf('Event %i out of %i',ii,nevent));
                for jj = 1:nrec
                    for kk = 1:ncomp
                        imatch = find( and( and( meas_array(:,2)==tt, meas_array(:,3)==ii), ...
                                            and( meas_array(:,5)==jj, meas_array(:,6)==kk) ));
                        if ~isempty(imatch)
                            % take the first only
                            wdiff_array(ii,jj,kk,tt,imod+1) = meas_array(imatch(1),16);
                        end
                    end
                end
            end
        end
    
    end
    
    save('wdiff_array.mat','wdiff_array');

else
   load('wdiff_array.mat'); 
end

whos wdiff_array

% total number of windows
N = length(meas_array);

%----------------------------------------------
% CODE IN PROGRESS
% for each record for each model (m0, m1) that has at least one window picked,
% we save the integrated waveform difference as the purest measure of misfit
% between the data and synthetics

[m1,m2,m3,m4,m5] = size(wdiff_array)

if 0==1

    ratio_array = []; ll = 0;
    for tt = 1:nper
        for ii = imin:imax
            disp(sprintf('Event %i out of %i',ii,nevent));
            for jj = 1:nrec
                for kk = 1:ncomp
                    % check that there is a measurement for all models
                    if prod(wdiff_array(ii,jj,kk,tt,:)) > 0
                        ll = ll+1;
                        rat = wdiff_array(ii,jj,kk,tt,imodel_max+1) / wdiff_array(ii,jj,kk,tt,imodel_min+1);
                        ratio_array(ll,:) = [tt ii jj kk rat];
                    end
                end
            end
        end
    end

    %ratio_sort = sortrows(ratio_array,5);
    ratio_sort = sortrows(ratio_array,[-1 5]);
    
    hmax = length(ratio_sort);
    hmax = 100;
    for hh = 1:hmax
        tt = ratio_sort(hh,1);
        ii = ratio_sort(hh,2);
        jj = ratio_sort(hh,3);
        kk = ratio_sort(hh,4);
       disp(sprintf('%10s --> %6.1f %7s %5s %6.3f',char(eids{ii}),...
           Tvec(tt),char(strec{jj}),char(comps{kk}),ratio_sort(hh,5) ));
    end
    
    %------------------------

    rec_array = []; ll = 0;
    for tt = 1:nper
        for ii = imin:imax
            disp(sprintf('Event %i out of %i',ii,nevent));
            for jj = 1:nrec

                    % check that there is a measurement for all components for all models
                    if prod( wdiff_array(ii,jj,:,tt,:) ) > 0
                        ll = ll+1;
                        rat = prod( wdiff_array(ii,jj,:,tt,imodel_max+1) ./ wdiff_array(ii,jj,:,tt,imodel_min+1) );
                        rec_array(ll,:) = [tt ii jj rat];
                    end
                end

        end
    end
    rec_sort = sortrows(rec_array,[-1 4]);
    
    for hh = 1:length(rec_sort)
        tt = rec_sort(hh,1);
        ii = rec_sort(hh,2);
        jj = rec_sort(hh,3);
       disp(sprintf('%10s --> %6.1f %7s %5s %10.6f',char(eids{ii}),...
           Tvec(tt),char(strec{jj}),rec_sort(hh,4) ));
    end
    
    %------------------------
    
%     % loop over all stations
%     tarray = []; ll = 0;
%     for jj = 1:nrec
%         imatch = find( ratio_array(:,3)==jj);
%         if length(imatch)==6
%             ll = ll+1;
%            tarray(ll,:) = [jj prod(ratio_array(imatch,5)) ];
%         end
%     end
%     tarray_sort = sortrows(tarray,2);
%     for hh = 1:20
%         jj = tarray_sort(hh,1);
%        disp(sprintf('%7s %6.3f',char(strec{jj}),tarray_sort(hh,2) ));
%     end

    figure; N = plot_histo(ratio_array(:,end),[0:0.2:5]);
    figure; N = plot_histo(log(ratio_array(:,end)),[-4:0.2:4]);

    break
end

%=======================================================================
