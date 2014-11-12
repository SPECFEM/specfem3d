%
% display_meas.m
% CARL TAPE, 13-Oct-2009
% printed xxx
%
% /cig/seismo/3D/mt_measure_adj_work/scripts_tomo/matlab/
%
% calls xxx
% called by xxx
%

function display_meas(meas_array,Tvec,eids,strec,comps)

X = meas_array;
[n,m] = size(X);
if m ~= 28, error(' should be 28 columns in the input matrix'); end

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

disp('----------------');
disp('INDEX (window, period, event, receiver, component)');
for ii = 1:n
   disp(sprintf('%10s --> %6.1f %7s %6s %3i DT = %6.2f +- %6.2f DA = %6.2f chi = %6.2f -- %i',...
       char(eids{index_event(ii)}),Tvec(index_per(ii)),char(strec{index_rec(ii)}),...
       char(comps{index_comp(ii)}),isub_win(ii),...
       measCC_dT(ii),sigmaCC_dT(ii),measCC_dA(ii),tr_chi(ii),index_win(ii) ));
   
   % optional: list MT measurment values
   %disp(sprintf('%40s DT = %6.2f +- %6.2f DA = %6.2f',' ',measMT_dT(ii),sigmaMT_dT(ii),measMT_dA(ii)));
end
disp('INDEX (window, period, event, receiver, component)');
disp('====================');

%======================================================
