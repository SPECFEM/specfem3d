#!/ichec/packages/octave/3.6.3/octave-3.6.3_build/bin/octave
clear all
close all
clc

sourceType = input('Please input the SOURCE type(FORWARD or BACKWARD) : ','s')

forwardSourceComponent=[0 0 1];
backwardSourceComponent=[1 1 1];

switch sourceType
case 'FORWARD'
[nt_status nt] = system('grep NSTEP ../in_data_files/Par_file | cut -d = -f 2');
[dt_status dt] = system('grep DT ../in_data_files/Par_file | cut -d = -f 2');
[f_status f] = system('grep half\ duration ../in_data_files/backup/CMTSOLUTION_FORWARD | cut -d : -f 2');

nt=str2num(nt);
dt=str2num(dt);
f=str2num(f);

sourceNumber = length(f);

s_z = zeros(nt,sourceNumber);
s_x = zeros(nt,sourceNumber);
s_y = zeros(nt,sourceNumber);

for nSource = 1:sourceNumber
  [tCut sCut] = ricker(f(nSource),dt);

  if forwardSourceComponent(1)
  s_x(1:length(sCut),nSource) = sCut;
  end
  if forwardSourceComponent(2)
  s_y(1:length(sCut),nSource) = sCut;
  end
  if forwardSourceComponent(3)
  s_z(1:length(sCut),nSource) = sCut;
  end
end


case 'BACKWARD'
[nt_status nt] = system('grep NSTEP ../in_data_files/Par_file | cut -d = -f 2');
nt=str2num(nt);

fid=fopen('../in_data_files/backup/STATIONS_FORWARD');
  c=textscan(fid,'%s %s %f %f %f %f');
fclose(fid);

stationNumber=length(c{1,1});

s_x=zeros(nt,stationNumber);
s_y=zeros(nt,stationNumber);
s_z=zeros(nt,stationNumber);

band='FX';
%component='Z';
variable='semd';

cutLength = 3000;
cutShift = 0;

totalCutTrace=zeros(cutLength,stationNumber);

benchStation=floor((1+stationNumber)/2);
bench=load(['../in_out_files/OUTPUT_FILES/' c{1,1}{benchStation} '.' c{1,2}{benchStation} '.' band 'Z' '.' variable]);
[benchMax benchMaxIndex] = max(bench(:,2));

for nStation = 1:stationNumber
trace_x = load(['../in_out_files/OUTPUT_FILES/' c{1,1}{nStation} '.' c{1,2}{nStation} '.' band 'X' '.' variable]);
trace_y = load(['../in_out_files/OUTPUT_FILES/' c{1,1}{nStation} '.' c{1,2}{nStation} '.' band 'Y' '.' variable]);
trace_z = load(['../in_out_files/OUTPUT_FILES/' c{1,1}{nStation} '.' c{1,2}{nStation} '.' band 'Z' '.' variable]);
cutTrace_x = trace_x(benchMaxIndex - cutLength/2+cutShift:benchMaxIndex + cutLength/2 -1+cutShift,2);
cutTrace_y = trace_y(benchMaxIndex - cutLength/2+cutShift:benchMaxIndex + cutLength/2 -1+cutShift,2);
cutTrace_z = trace_z(benchMaxIndex - cutLength/2+cutShift:benchMaxIndex + cutLength/2 -1+cutShift,2);
totalCutTrace_x(:,nStation) = cutTrace_x;
totalCutTrace_y(:,nStation) = cutTrace_y;
totalCutTrace_z(:,nStation) = cutTrace_z;
end

timeReversedTotalCutTrace_x = flipud(totalCutTrace_x);
timeReversedTotalCutTrace_y = flipud(totalCutTrace_y);
timeReversedTotalCutTrace_z = flipud(totalCutTrace_z);
figure
plot(timeReversedTotalCutTrace_x)
figure
plot(timeReversedTotalCutTrace_y)
figure
plot(timeReversedTotalCutTrace_z)
pause(20)

if backwardSourceComponent(1)
s_x(1:cutLength,1:stationNumber) = timeReversedTotalCutTrace_x;
end
if backwardSourceComponent(2)
s_y(1:cutLength,1:stationNumber) = timeReversedTotalCutTrace_y;
end
if backwardSourceComponent(3)
s_z(1:cutLength,1:stationNumber) = timeReversedTotalCutTrace_z;
end

otherwise
error('Wrong SOURCE type.')
end
save("-ascii",['../in_data_files/backup/SOURCE_X_' sourceType],'s_x')
save("-ascii",['../in_data_files/backup/SOURCE_Y_' sourceType],'s_y')
save("-ascii",['../in_data_files/backup/SOURCE_Z_' sourceType],'s_z')
