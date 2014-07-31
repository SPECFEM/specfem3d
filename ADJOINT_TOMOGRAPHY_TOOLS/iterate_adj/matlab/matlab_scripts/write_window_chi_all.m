%
% write_window_chi_all.m
% CARL TAPE, 30-July-2009
% printed xxx
%
% This reads in an array of window_chi misfit files and outputs a text file
% that can be used for additional analysis and plots (e.g., GMT histograms).
%
% calls xxx
% called by xxx
%

function write_window_chi_all(X,ofile,stmod,Tminvec,Tmaxvec,eids,strec,stnet,comps)

% total number of measurements
n = length(X);

index_win    = X(:,1);
index_per    = X(:,2);
index_event  = X(:,3);
index_net    = X(:,4);
index_rec    = X(:,5);
index_comp   = X(:,6);
isub_win     = X(:,7);
iker         = X(:,8);
%---------------------------
seisdur      = X(:,9);
windur       = X(:,10);
T_pmax_dat   = X(:,11);
T_pmax_syn   = X(:,12);
%---------------------------
seis_d2      = X(:,13);
win_d2       = X(:,14);
seis_s2      = X(:,15);
win_s2       = X(:,16);
seis_diff2   = X(:,17);
win_diff2    = X(:,18);
%---------------------------
measCC_dT    = X(:,19);
sigmaCC_dT   = X(:,20);
measCC_dA    = X(:,21);
sigmaCC_dA   = X(:,22);
measMT_dT    = X(:,23);
sigmaMT_dT   = X(:,24);
measMT_dA    = X(:,25);
sigmaMT_dA   = X(:,26);
tr_chi       = X(:,27);
am_chi       = X(:,28);

%n = 1000;    % testing

fid = fopen(ofile,'w');

fprintf(fid,['%4s%10s%7s%7s%9s%5s%3s%7s%7s%6s%6s' ...
        repmat('%10s',1,6) repmat('%7s',1,10) '\n'],...
        'mod','eid','Tmin','Tmax','rec','comp','iw',...
        'sdur','wdur','Tdat','Tsyn','sd2','wd2','ss2','ws2','sd2','wd2',...
        'DTCC','sigma','DACC','sigma','DTMT','sigma','DAMT','sigma',...
        'tr_chi','am_chi');

for ii = 1:n
    fprintf(fid,['%4s%10s%7.1f%7.1f%9s%5s%3i%7.1f%7.1f%6.1f%6.1f' ...
        repmat('%10.2e',1,6) repmat('%7.3f',1,10) '\n'],...
       stmod,char(eids{index_event(ii)}),Tminvec(index_per(ii)),Tmaxvec(index_per(ii)),...
       char(strec{index_rec(ii)}),char(comps{index_comp(ii)}),isub_win(ii),...
       X(ii,9:12),X(ii,13:18),X(ii,19:28) );
end
fclose(fid);

%------------

% find MT measurements
nMT = 0;      % temporary

% values to go with the histograms
%  1 - nwinCC
%  2 - m16 DT-CC mean
%  3 - m16 DT-CC std
%  4 - m16 DT-MT mean
%  5 - m16 DT-MT std
%  6 - nwinMT
%  7 - m16 DlnA mean
%  8 - m16 DlnA std
%  9 - m16 DlnA-MT mean
% 10 - m16 DlnA-MT std

disp('writing out statistic values...');
fid = fopen([ofile '_stats'],'w');
fprintf(fid,'%10i%12.4e%12.4e%12.4e%12.4e%10i%12.4e%12.4e%12.4e%12.4e',...
    n,mean(measCC_dT),std(measCC_dT),0,0,nMT,mean(measCC_dA),std(measCC_dA),0,0);
fclose(fid);

%======================================================
