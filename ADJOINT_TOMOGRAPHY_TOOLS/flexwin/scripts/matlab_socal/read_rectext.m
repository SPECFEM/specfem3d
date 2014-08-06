%
% function [zrt_win,name,netwk] = read_rectext(filename)
% CARL TAPE, 27-July-2007
% printed xxx
%
% This reads in a particular file output from the windowing code, and it
% returns the number of windows picked per trace.
%
% calls xxx
% called by window_summary.m
%

function [zrt_win,name,netwk] = read_rectext(filename)

[x,y,junk1,junk2,junk3,junk4,text1,text2,text3,text4,text5,suffix] ...
   = textread(filename,'%f%f%f%f%f%s%s%s%s%s%s%s','headerlines',0);

dist = y;           % arc-distance
name = text1;       % station name
netwk = text2;      % station network
n = length(x);      

zrt_win = zeros(n,3);
for ii=1:n
    st1 = char(text3{ii});
    st2 = char(text4{ii});
    st3 = char(text5{ii});
    stcat{ii} = [st1 st2 st3];

    zrt_win(ii,:) = [str2num(st1(2)) str2num(st2(1)) str2num(st3(1))];
end
%win_tot_vec = sum(zrt_win,2);  % total window picks for each receiver
%win_tot = sum(win_tot_vec);

%======================================================