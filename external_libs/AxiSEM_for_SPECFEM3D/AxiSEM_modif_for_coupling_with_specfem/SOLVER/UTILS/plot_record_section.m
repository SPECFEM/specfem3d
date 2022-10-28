clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load info_matlab.dat
nrec=info_matlab(1); nt=info_matlab(2); dt=info_matlab(3); t0=info_matlab(4);
ii=5;
for i=1:nrec; epidist(i)=info_matlab(ii)*180./pi;  ii=ii+1; end
clear info_matlab
filen=sprintf('info_seiscomp.dat');fid=fopen(filen,'r'); reccomp=fscanf(fid,'%s ',[3]); reccomp=reccomp'; %'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   LOAD SEISMOGRAMS   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filen=sprintf('info_seisnames.dat'); fid=fopen(filen,'r');
cd('SEISMOGRAMS')
for j=1:3;  for i=1:nrec
      seisfile=fgets(fid,100); filen=sprintf(strtrim(seisfile)); seistmp=load(strtrim(filen)); seis(:,i,j)=seistmp(:,2);
      if j==1 && i==1; time=seistmp(:,1); end
end; end
clear seistmp
cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   PLOT TRAVELTIMES (TAUP)   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
tau_exist=exist('TAUP/taup_P_traveltime.dat','file');
if tau_exist > 0
load('TAUP/taup_PcP_traveltime.dat'); load('TAUP/taup_pP_traveltime.dat');load('TAUP/taup_S_traveltime.dat')
load('TAUP/taup_PP_traveltime.dat'); load('TAUP/taup_P_traveltime.dat')
load('TAUP/taup_ScS_traveltime.dat'); load('TAUP/taup_SS_traveltime.dat')
if max(epidist) >95. ;load('TAUP/taup_Sdiff_traveltime.dat'); load('TAUP/taup_Pdiff_traveltime.dat'); load('TAUP/taup_SKS_traveltime.dat'); end
if max(epidist) >95. ; load('TAUP/taup_SKKS_traveltime.dat'); end
for j=1:3
   figure(j)
   plot(taup_P_traveltime(:,2),taup_P_traveltime(:,1),'r');  hold on
   plot(taup_pP_traveltime(:,2),taup_pP_traveltime(:,1),'b'); plot(taup_PP_traveltime(:,2),taup_PP_traveltime(:,1),'r--')
   plot(taup_PcP_traveltime(:,2),taup_PcP_traveltime(:,1),'m-.'); plot(taup_S_traveltime(:,2),taup_S_traveltime(:,1),'c')
   plot(taup_SS_traveltime(:,2),taup_SS_traveltime(:,1),'b--'); plot(taup_ScS_traveltime(:,2),taup_ScS_traveltime(:,1),'g-.')
   plot(taup_SKS_traveltime(:,2),taup_SKS_traveltime(:,1),'g--'); plot(taup_SKKS_traveltime(:,2),taup_SKKS_traveltime(:,1),'b--')
   legend('P','pP','PP','PcP','S','SS','ScS','SKS')
   if max(epidist) >100.
     plot(taup_Sdiff_traveltime(:,2),taup_Sdiff_traveltime(:,1),'c-.'); plot(taup_Pdiff_traveltime(:,2),taup_Pdiff_traveltime(:,1),'r-.')
end; end; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   PLOT SEISMOGRAMS   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'DefaulttextInterpreter', 'none')
deltadist=max(epidist)-min(epidist);
filen2=sprintf('info_seisstations.dat'); fid2=fopen(filen2,'r');
for j=1:3
  globmax(j)=max(max(abs(seis(:,:,j)))); if globmax(j)==0. ; globmax(j)=1.; end
  figure(j); for i=1:nrec
       plot(seis(:,i,j)/globmax(1)*epidist(i)^1.5/(epidist(2))*deltadist/nrec+epidist(i),time,'k')
       hold on
      % plot station names
      textfile=fgets(fid2,100); filen2=sprintf(strtrim(textfile)); text(epidist(i),time(1),sprintf('%3s',filen2))
   end
   ylabel('time [s]'); xlabel('epicentral distance [^o]'); title(sprintf('displacement component: %s',reccomp(j)))
   axis([min(epidist)-3*deltadist/nrec max(epidist)+3*deltadist/nrec min(time) max(time) ])
end
