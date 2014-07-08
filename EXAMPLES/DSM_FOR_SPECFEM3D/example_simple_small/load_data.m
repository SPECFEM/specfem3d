clc
clear all
format long

N = input('please input the total number of the traces for plot:');
%if (N <10)
%N= 10;
%end 
   
element1 = 'X'
element2 = 'Y'
element3 = 'Z'

name = 'v'
S = 1;

tic

index2 = strcat('S0000',num2str((1:9)','%-d'),'.TS','.BX',num2str(element1),'.sem',name);
index3 = strcat('S0000',num2str((1:9)','%-d'),'.TS','.BX',num2str(element2),'.sem',name);
index4 = strcat('S0000',num2str((1:9)','%-d'),'.TS','.BX',num2str(element3),'.sem',name);
if (N > 9)
%   A = min (N-10,9);
index21 = strcat('S000',num2str((10:N)','%-d'),'.TS','.BX',num2str(element1),'.sem',name);
index31 = strcat('S000',num2str((10:N)','%-d'),'.TS','.BX',num2str(element2),'.sem',name);
index41 = strcat('S000',num2str((10:N)','%-d'),'.TS','.BX',num2str(element3),'.sem',name);
end
 
if (N >99)
N=99;
end


order = size(index2);

if (N>9)
index2(10:N,1:order(2))=0;
index2(10:N,:)=index21;

index3(10:N,1:order(2))=0;
index3(10:N,:)=index31;

index4(10:N,1:order(2))=0;
index4(10:N,:)=index41;
end

clear index21 index31 index41 order


mid = load(index2(S,:));
order = size(mid);

data = zeros(order(1),order(2),N);

for k=S:1:N
    data(:,:,k)= load(index2(k,:));
end

y2 = zeros(order(1),N);
for k=S:1:N
    y2(:,k) = data(:,2,k);
end

for k=S:1:N
    data(:,:,k)= load(index3(k,:));
end

y3 = zeros(order(1),N);
for k=S:1:N
    y3(:,k) = data(:,2,k);
end

for k=S:1:N
    data(:,:,k)= load(index4(k,:));
end

y4 = zeros(order(1),N);
for k=S:1:N
    y4(:,k) = data(:,2,k);
end
%y2 = y2*1d6;
y3 = y3*1d6;
%y4 = y4*1d6;

clear mid 

trace_s = input('please input the start trace for plot:');
trace_e = input('please input the end trace for plot:');
if (trace_s < S)
    trace_s = S;
end

if (trace_s < 0 | trace_e < 0)
trace_s = S
trace_e = N
elseif (trace_e > N | trace_s > N)
 if (trace_s >N )
 trace_s = S
 elseif (trace_s < S)
 trace_s = S
 else
 trace_s = floor(trace_s)
 end
trace_e = N
else
 if (trace_s < S)
 trace_s = S
 else
 trace_s = floor(trace_s)
 end
trace_e = floor(trace_e)
end


if (trace_s > trace_e)
    middle = trace_s;
    trace_s=trace_e;
    trace_e=middle;
end


max_value(1) = max(max(abs(y2(:,trace_s:trace_e))));
max_value(2) = max(max(abs(y3(:,trace_s:trace_e))));
max_value(3) = max(max(abs(y4(:,trace_s:trace_e))))
max_amp = max(abs(y3(:,trace_s:trace_e)));
max_amp_xz = max(abs(y3(:,trace_s:trace_e)));
disp = max ( max_value(:) )/2.0d0
%disp = 0.0d0

'error between 2 results is:'
%max(max(abs(y2(:,trace_s:trace_e)-y3(:,trace_s:trace_e))))/max(max1+max2)
%disp=0;


scale = input('please input the scaling factor which is a integer:');
%scale = 3; !!scale is the scaling ratio for the plot data compression
col = floor(order(1)/scale);
x=zeros(col,N);
y=x;
z=x;
for k=S:1:N
    x(:,k)=y2(scale:scale:order(1),k);
    y(:,k)=y3(scale:scale:order(1),k);
    z(:,k)=y4(scale:scale:order(1),k);
end

clear y2 y3 y4

xxx=data(scale:scale:order(1),1,S);
%xxx=xxx+10.0d0;
order1 = size(data);
clear data

ts = input('please input the start time for plot:');
te = input('please input the end time for plot:');

if (ts < 0 | te <0)
order = size(x);
ts = 1;
te = order(1);
elseif (te > order1(1)*0.1)
ts = floor(ts/0.1/scale) +1; 
order = size(x);
te = order(1);
else
ts = floor(ts/0.1/scale) +1;
te = floor(te/0.1/scale) -1;
end


zoomin = input('Please input the magnifying factor:');
if (zoomin <= 0)
    zoomin = 1;
else
    x=x*zoomin;
    y=y*zoomin;
    z=z*zoomin;
end

error = y-z;
error_xz = x-z;
disp_error =  max( max(abs(y(:,trace_s:trace_e)-z(:,trace_s:trace_e))) )
disp_error_xz = max( max(abs(x(:,trace_s:trace_e)-z(:,trace_s:trace_e))) )
sumd = zeros(te-ts+1,trace_e);
angular = zeros(te-ts+1,trace_e);
for k=trace_s:1:trace_e
    sumd(:,k) = ( x(ts:te,k).^2 + y(ts:te,k).^2 + z(ts:te,k).^2 ).^(0.5);
    angular(:,k) = acos( y(ts:te,k)./sumd(1:te-ts+1,k) );
end


for k=trace_s:1:trace_e
    x(:,k) = x(:,k) + (k-1)*disp;
    y(:,k) = y(:,k) + (k-1)*disp;
    z(:,k) = z(:,k) + (k-1)*disp;
    error(:,k) = error(:,k) + (k-1)*disp_error;
    error_xz(:,k) = error_xz(:,k) + (k-1)*disp_error_xz;
end 

'output the reference relative error of trace 2 and 3:'
%disp_error/max ( max_value(:) )
 error_array = max(abs(y(:,trace_s:trace_e)-z(:,trace_s:trace_e))) ./ max_amp 
 relative_error = max(error_array)
 error_array_xz = max(abs(x(:,trace_s:trace_e)-y(:,trace_s:trace_e))) ./ max_amp_xz 
 relative_error_xz = max(error_array_xz)


option = input('please input the option for plot:(1 for element1, 2 for element2, 3 for element 3, 4 for absolute value, 0 for error between Y and Z trace, positive number for all, negtive for exit)  ');

while(option >= 0)
if (option ==1) 
    plot(x(ts:te,trace_s:trace_e),xxx(ts:te),'bd-.','linewidth',1,'MarkerSize',3);
elseif (option ==2)
    plot(y(ts:te,trace_s:trace_e),xxx(ts:te),'rs:','linewidth',1,'MarkerSize',3);
elseif (option ==3)
    plot(z(ts:te,trace_s:trace_e),xxx(ts:te),'gp:','linewidth',1,'MarkerSize',3);
    
elseif (option ==4)
    plot(sumd(1:te-ts+1,trace_s:trace_e),xxx(ts:te),'c+:','linewidth',1,'MarkerSize',3);
elseif (option ==0) %error of the 2 data sets
    plot(error(ts:te,trace_s:trace_e),xxx(ts:te),'kd:','linewidth',1,'MarkerSize',2);
elseif (option ==23) %comparation of the 2 data sets
    plot(y(ts:te,trace_s:trace_e),xxx(ts:te),'rs:',z(ts:te,trace_s:trace_e),xxx(ts:te),'gp:','linewidth',1,'MarkerSize',3);    
elseif (option ==12) %comparation of the 2 data sets
    plot(x(ts:te,trace_s:trace_e),xxx(ts:te),'bd-.',y(ts:te,trace_s:trace_e),xxx(ts:te),'rs:','linewidth',1,'MarkerSize',3);    
elseif (option ==13) %comparation of the 2 data sets
    plot(x(ts:te,trace_s:trace_e),xxx(ts:te),'bd-.',z(ts:te,trace_s:trace_e),xxx(ts:te),'gp:','linewidth',1,'MarkerSize',3);    

else
    plot(x(ts:te,trace_s:trace_e),xxx(ts:te),'bd-.',y(ts:te,trace_s:trace_e),xxx(ts:te),'rs:',z(ts:te,trace_s:trace_e),xxx(ts:te),'gp:','linewidth',1,'MarkerSize',3);
end
title('data of wavefields');
ylabel('Time (s)');
xlabel('displacement blue for element1 / red for elemnt2 / green for element3');
option = input('please input the option for plot:(1 for element1, 2 for element2, 3 for element 3, 4 for absolute value, positive number for all, negtive for exit)  ');
end

%  matlab plot color and type of the point  
%     y         yellow           ·              
%     m         pink             ○              
%     c         bright blue      ×              
%     r         red              +           
%     g         green            -             
%     b         blue             *                
%     w         white            :           
%     k         black            -.
                             

toc
