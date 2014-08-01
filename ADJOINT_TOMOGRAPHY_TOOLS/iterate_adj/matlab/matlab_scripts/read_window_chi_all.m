%
% read_window_chi_all.m
% Carl Tape, 14-March-2008
%
% This reads in a set of measurement files that are produced by the
% measure_adj softward package for making measurements and adjoint sources.
%
% calls read_window_chi.m
% called by compute_misfit.m, compare_misfit.m
%

function meas_array = read_window_chi_all(imin,imax,stmod,Tminvec,Tmaxvec,dir0,eids,strec,stnet,comps)

nper = length(Tminvec);
nevent = length(eids);
nnet = length(stnet);
nrec = length(strec);
ncomp = length(comps);

if imax > nevent, error('imax is greater than nevent'); end

k2 = 0;
meas_array = [];
%for tt = 1:1
for tt = 1:nper
    Tmin = Tminvec(tt);
    Tmax = Tmaxvec(tt);
    Ttag = sprintf('T%3.3i_T%3.3i',Tmin,Tmax);
    Ttag2 = [Ttag '_' stmod];
    disp('-----------------------------------------');

    for ii = imin:imax    % loop over events

        disp(sprintf('Event %i out of %i -- %s',ii,nevent,eids{ii}));
        dir1 = [dir0 eids{ii} '/' stmod '/MEASURE_' Ttag '/'];
        mfile = [eids{ii} '_' Ttag2 '_window_chi'];
        filename = [dir1 mfile];

        if ~exist(filename,'file')
            disp('--> file does not exist (or nwin = 0)');
        else
            % read the window_chi file
            [stnet0,strec0,comp,iwin,iker,t1,t2,...
            chiMT_dT,chiMT_dA,chiCC_dT,chiCC_dA,...
            measMT_dT,measMT_dA,measCC_dT,measCC_dA,...
            sigmaMT_dT,sigmaMT_dA,sigmaCC_dT,sigmaCC_dA,...
            wind2,wins2,windiff2,windur, ...
            seisd2,seiss2,seisdiff2,seisdur,...
            tr_chi,am_chi,T_pmax_dat,T_pmax_syn] = read_window_chi(filename);

            % waveform differences normalized by data-squared
            %win_chi = windiff2 ./ wind2;        % window only
            %seis_chi = seisdiff2 ./ seisd2;        % full record (repeated for seismograms with more than one window)

            % a waveform measure of the power in a window relative to the full record
            %win_reld2 = wind2 ./ seisd2;
            %win_rels2 = wins2 ./ seiss2;

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

            % NOTE: remove the factors of 0.5 that is used in the wdiff
            %       definition in measure_adj.f90
            meas_array(kinds,:) = [kinds tt*ones(nwin,1) ii*ones(nwin,1) index_net index_rec index_comp iwin iker ...
                seisdur windur T_pmax_dat T_pmax_syn ...
                2*seisd2  2*wind2 2*seiss2  2*wins2 2*seisdiff2 2*windiff2 ...
                measCC_dT sigmaCC_dT measCC_dA sigmaCC_dA ...
                measMT_dT sigmaMT_dT measMT_dA sigmaMT_dA tr_chi am_chi];
            %  1  kinds
            %  2  index_T
            %  3  index_event
            %  4  index_network
            %  5  index_rec
            %  6  index_comp
            %  7  iwin
            %  8  iker
            %  9  seisdur
            % 10  windur
            % 11  T_pmax_dat
            % 12  T_pmax_syn
            %-------------------
            % waveform difference information
            % 13  seisd2
            % 14  wind2
            % 15  seiss2
            % 16  wins2
            % 17  seisdiff2
            % 18  windiff2   
            %-------------------
            % other misfit measures
            % 19  measCC_dT
            % 20  sigmaCC_dT
            % 21  measCC_dA
            % 22  sigmaCC_dA
            % 23  measMT_dT
            % 24  sigmaMT_dT
            % 25  measMT_dA
            % 26  sigmaMT_dA
            % 27  tr_chi
            % 28  am_chi
        end

    end
end

%==========================================================================
