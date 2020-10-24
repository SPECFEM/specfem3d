function [ ] = NOISE_TOMOGRAPHY(NSTEP,dt,Tmin,Tmax,NOISE_MODEL)
% This is a short Matlab program used to sample the spectrum of the 
%     Peterson's noise model, and convert it to time-domain --- the 
%     source time function we will use in NOISE TOMOGRAPHY simulations.
%
%***********************************************************************
%*******Please read the manual of SPECFEM3D package for guidance********
%***********************************************************************
%
% Usage:
%   NOISE_TOMOGRAPHY(NSTEP,dt,Tmin,Tmax,NOISE_MODEL)
%
% Example:
%   NOISE_TOMOGRAPHY(2999,0.05,10.,20.,'NLNM')
%
% /////////////////////////////////////////////////////////////////////////
% Parameters:
% NSTEP       --- number of time steps (always odd for NOISE TOMOGRAPHY)
% dt          --- time interval of specfem3D solver
% Tmin, Tmax  --- the period range you are working with (in seconds)
% NOISE_MODEL --- the Peterson's noise model, 
%                 either 'NLNM' or 'NHNM' (with the quote!)
%                 'NLNM': New Low  Noise Model (in 1993, the model was New)
%                 'NHNM': New High Noise Model
%                 or 'FLAT': for a flat noise spectrum (uniform noise within period range)
%
% /////////////////////////////////////////////////////////////////////////
%
% ATTENTION:
% ***NEVER*** try to calculate "NSTEP" and "dt" by yourself!
% They can be found when you compile the SPECFEM3D package,
% with the correct DATA/Par_file
% If DATA/Par_file is not specified for noise tomography simulations,
% you won't get correct values either!
% 
% It is highly recommended that you compile the package with the correct 
% DATA/Par_file you will be using
%
% 
% Output from the example above:
% a figure + following message 
% (the path /data2/yangl/3D_NOSIE/ must be different)
% *************************************************************
% the source time function has been saved in:
% /data2/yangl/3D_NOISE/S_squared
% S_squared should be put into directory:
% ./NOISE_TOMOGRAPHY/ in the SPECFEM3D package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER PARAMETERS
%% figure plotting
show_figures = true;

%% taper option
taper_type = 1;                % cosine type (1=on/0==off)
taper_length_percentage = 0.1; % taper length as a percentage of full length
taper_length_min = 40;         % minimum number of steps for tapering ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if show_figures
  clf;
  fontsize=10;       % fontsize for figures
end

%% derived parameters
T=(NSTEP-1)*dt;    % total simulation time
N_mid=floor((NSTEP+1)/2); % (NSTEP+1)/2 is the middle of the (NSTEP-1) points
fmax=1/2/dt;       % Nyquist frequency, due to sampling theory
df=1/T;            % frequency interval
f=[0:df:fmax -fmax:df:-df]; % discrete frequencies

%% checks length
if T<Tmax
    fprintf(['The simulation length T(' num2str(T) ...
        ') is smaller than the required maximum period Tmax(' ...
        num2str(Tmax) ')\n']);
    return;
end;

%% checks noise model string
use_flat_noise_spectrum = false;
if NOISE_MODEL=='NLNM'
  model_info='Peterson New Low Noise Model';
elseif NOISE_MODEL=='NHNM'
  model_info='Peterson New High Noise Model';
elseif NOISE_MODEL=='FLAT'
  % use a flat noise spectrum instead of Peterson model
  use_flat_noise_spectrum = true;
  model_info='uniform noise spectrum';
else
  fprintf('Error: noise model %s not recognized, use NLNM, NHNM or FLAT for low, high or flat noise model',NOISE_MODEL);
  exit 1;
end

%% user output
fprintf('NOISE_TOMOGRAPHY input:\n');
fprintf('  number of time steps = %i \n',NSTEP);
fprintf('  time step size       = %f s\n',dt);
fprintf('  period range min/max = %f / %f s\n',Tmin,Tmax);
fprintf('  noise model          = %s - %s\n',NOISE_MODEL,model_info);
fprintf('\n');
fprintf('  total simulation time = %f s\n',T);
fprintf('  Nyquist frequency: %f Hz\n',fmax);

%% initialize the power spectrum of noise
accel=zeros(size(f)); % for acceleration
veloc=zeros(size(f)); % for velocity
displ=zeros(size(f)); % for displacement
filter_array=zeros(size(f)); % for constraining frequencies within ranges
                             % defined by Tmin & Tmax

%% calculate the power spectrum of noise from Peterson's model (1993)
fprintf('\n  calculating noise power spectrum: number of positive frequencies = %i \n',N_mid);
% only calculate for positive frequencies
% the negative frequencies will be updated later using symmetry
for l=1:N_mid
  if use_flat_noise_spectrum
     accel(l)=1; veloc(l)=1; displ(l)=1;
  else
    [accel(l) veloc(l) displ(l)]=PetersonNoiseModel(1/f(l),NOISE_MODEL);
  end
end

% plots figure
if show_figures
  fprintf('  plotting figure: Power spectrum\n');
  figure(1);
  subplot(2,2,1);
  semilogx(f(1:N_mid),accel(1:N_mid),'b');hold on;
  semilogx(f(1:N_mid),veloc(1:N_mid),'g');hold on;
  semilogx(f(1:N_mid),displ(1:N_mid),'r');hold on;
  legend('acceleration','velocity','displacement');
  xlabel('Frequency (Hz)','fontsize',fontsize);ylabel('Amplitude (dB)','fontsize',fontsize);
  title('Power Spectrum of Peterson''s Noise Model in dB','fontsize',fontsize);
  xlim([1e-4 1e1]);ylim([-250 -50]);
  %xlim([0.9*df 1.1*fmax]); ylim([-250 -50]);
end

%% change power spectrum from dB to real physical unit
for l=1:N_mid
        accel(l)=10.^(accel(l)/10);
        veloc(l)=10.^(veloc(l)/10);
        displ(l)=10.^(displ(l)/10);
end

%debug
%n = size(displ,2);tmp=zeros(n,2);
%tmp(:,1)=[1:n];tmp(:,2)=displ(:);
%save('tmp_displ.txt','tmp','-ASCII');

%% constrain the power spectrum only within the range [Tmin Tmax]
fprintf('\n  filtering power spectrum:\n    period     min/max = %f / %f\n    frequency  min/max = %f / %f \n',Tmin,Tmax,1/Tmax,1/Tmin);
for l=1:N_mid
    if abs(f(l))>=1/Tmax && abs(f(l))<=1/Tmin
      if use_flat_noise_spectrum
        filter_array(l)=1;
      else
        filter_array(l)=sin((f(l)-1/Tmax)/(1/Tmin-1/Tmax)*pi);
      end
    elseif abs(f(l))>=1/1.5/Tmax && abs(f(l))<1/Tmax
        filter_array(l)=sin ((f(l)-1/1.5/Tmax)/(1/Tmax-1/1.5/Tmax)*pi/2);
    elseif abs(f(l))<=1.5/Tmin && abs(f(l))>1/Tmin
        filter_array(l)=sin (-(f(l)-1/Tmin)/(1.5/Tmin-1/Tmin)*pi/2)+1;
    end
end

%debug
%n=size(filter_array,2); tmp=zeros(n,2);
%tmp(:,1)=[1:n];tmp(:,2)=filter_array(:);
%save('tmp_filter.txt','tmp','-ASCII');


accel=accel.*filter_array;
veloc=veloc.*filter_array;
displ=displ.*filter_array;

% plots figure
if show_figures
  fprintf('  plotting figure: Filtered Power spectrum\n');
  subplot(2,2,2);
  plot(f(1:N_mid),accel(1:N_mid)/max(abs(accel)),'b');hold on;
  plot(f(1:N_mid),veloc(1:N_mid)/max(abs(veloc)),'g');hold on;
  plot(f(1:N_mid),displ(1:N_mid)/max(abs(displ)),'r');hold on;
  xlabel('Frequency (Hz)','fontsize',fontsize);ylabel('Amplitude','fontsize',fontsize);
  title(['Power Spectrum filtered between [' num2str(Tmin) ' ' num2str(Tmax) '] s'],'fontsize',fontsize);
  legend(['accleration, scaled by ', num2str(max(abs(accel))), ' m^2/s^4/Hz'], ...
         ['velocity, scaled by ', num2str(max(abs(veloc))), ' m^2/s^2/Hz'],     ...
         ['displacement, scaled by ', num2str(max(abs(displ))), ' m^2/Hz']);
  xlim([0.8/Tmax 1.2/Tmin]);ylim([-0.1 1.5]);
  %xlim([0 min(1.5*Tmax,1.2/df)]);ylim([-0.1 1.5]);
end

% using only displacement from here on
clear accel;
clear veloc;

%% update power spectrum in the negative frequencies, using symmetry
% note the power spectrum is always REAL, instead of COMPLEX
for l=N_mid+1:NSTEP
    %accel(l)=conj(accel(NSTEP-l+2));
    %veloc(l)=conj(veloc(NSTEP-l+2));
    displ(l)=conj(displ(NSTEP-l+2));
end

%debug
%n=size(displ,2); tmp=zeros(n,2);
%tmp(:,1)=[1:n];tmp(:,2)=displ(:);
%save('tmp_displ2.txt','tmp','-ASCII');


%% prepare source time function for ensemble forward source -- S_squared
fprintf('\n  preparing source time function S_squared:\n    NSTEP = %i / dt = %f \n',NSTEP,dt);

% the file S_squared should be put into directory ./NOISE_TOMOGRAPHY/
% together with other two files: irec_main_noise & nu_main
S_squared=zeros(NSTEP,2); % first column: time (not used in SPECFEM3D package)
                          % second column: source time function
S_squared(:,1)=([1:NSTEP]-N_mid)*dt;
S_squared(:,2)=real(ifft(displ));

%debug
%save('tmp_S.txt','S_squared','-ASCII');

% change the order of the time series
% instead of having t=[0 T], we need t=[-T/2 T/2];
temp=S_squared(:,2);
temp(N_mid:NSTEP)=S_squared(1:N_mid,2);
temp(1:N_mid-1)  =S_squared(N_mid+1:NSTEP,2);

%debug
%n=size(temp,1); tmp=zeros(n,2);
%tmp(:,1)=[1:n];tmp(:,2)=temp(:);
%save('tmp_temp1.txt','tmp','-ASCII');  
  
% tapers ends
if taper_type==1
  fprintf('\n  tapering source time function S_squared\n');

  % determines taper length
  taper_length = round(taper_length_percentage * NSTEP);
  if taper_length < taper_length_min
    taper_length = taper_length_min;
  end
  fprintf('  using taper (cosine) length: %i \n',taper_length);
  % check size
  if NSTEP <= 2*taper_length
    fprintf('\nError number of time steps = %i must be bigger than taper window size = %i ! Please retry... \n',NSTEP,2*taper_length+1);
    exit(1);
  end

  % cosine taper (both branches, value 0 at index 1, value 1 at index length+1, value 0 at 2*length+1
  taper=zeros(2*taper_length+1);
  for l=1:2*taper_length+1
    taper(l) = (1.0-cos(pi*2.0*(l-1)/(2*taper_length)))/2.0;
  end

  % sets up window function for multiplication with signal
  % note: this becomes too memory intensive for large number of time steps, using taper directly
  %window=ones(NSTEP);
  %% adds increasing branch
  %window(1:taper_length) = taper(1:taper_length);
  %% adds decreasing branch
  %window(NSTEP-taper_length+1:NSTEP) = taper(taper_length+2:2*taper_length+1);
  %  
  %for l=1:NSTEP
  %  temp(l) = temp(l) * window(l);
  %end

  %debug
  %n=size(taper,2); tmp=zeros(n,2);
  %tmp(:,1)=[1:n];tmp(:,2)=taper(:,1);
  %save('tmp_taper.txt','tmp','-ASCII');

  % applies taper
  for l=1:taper_length
    % increasing branch
    temp(l) = temp(l) * taper(l);
    % decreasing branch
    k = NSTEP-taper_length+l;
    temp(k) = temp(k) * taper(taper_length+1+l);
  end

  %debug
  %n=size(temp,1); tmp=zeros(n,2);
  %tmp(:,1)=[1:n];tmp(:,2)=temp(:);
  %save('tmp_temp2.txt','tmp','-ASCII');  

  % plots figure
  if show_figures
    fprintf('  plotting figure: Taper and Taper window \n');
    figure(2);
    % taper
    subplot(2,1,1);
    plot(taper);hold on;
    xlabel('steps','fontsize',fontsize);ylabel('Taper','fontsize',fontsize);
    title('Taper','fontsize',fontsize);
    % taper window
    %subplot(2,1,2);
    %plot(window);
    %xlabel('steps','fontsize',fontsize);ylabel('Taper window','fontsize',fontsize);
    %title('Window','fontsize',fontsize);
    drawnow;
  end
end

% stores source time function values
S_squared(:,2)=temp;

% % filter the source time function if needed
% Wn=[1/Tmax 1/Tmin]/fmax;
% [B,A] = butter(4,Wn);
% S_squared(:,2)=filter(B,A,S_squared(:,2));

% plots figure
if show_figures
  fprintf('  plotting figure: Source time function S_squared \n');
  subplot(2,2,3);
  plot(S_squared(:,1)/60,S_squared(:,2),'r');
  xlabel('Time (min)','fontsize',fontsize); ylabel('Amplitude','fontsize',fontsize);
  xlim([-T/2 T/2]/60);
  title('Source Time Function for Ensemble Forward Source','fontsize',fontsize);

  subplot(2,2,4);
  plot(S_squared(:,1)/60,S_squared(:,2),'r');
  xlabel('Time (min)','fontsize',fontsize); ylabel('Amplitude','fontsize',fontsize);
  xlim([-Tmax Tmax]*0.10/60);
  %xlim([-1 1]/60*min(T/2,5*Tmax));
  title('Zoom-in of the Source Time Function','fontsize',fontsize);
  drawnow;
end

%% output the source time function
fprintf('\n  saving S_squared as ASCII file\n');
filename = 'S_squared';
save(filename,'S_squared','-ASCII')

%% user output
DIR=pwd;
fprintf('\n*************************************************************\n');
fprintf('the source time function has been saved in:\n');
fprintf([DIR '/S_squared\n']);
fprintf('\nS_squared should be put into directory:\n');
fprintf('./NOISE_TOMOGRAPHY/ in the SPECFEM3D package\n');
fprintf('*************************************************************\n');

% stop execution
%quit force;


