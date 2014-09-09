%
% readSTATIONS.m
% CARL TAPE, 15-March-2007
% printed xxx
%
% This function inputs a STATIONS file formatted for SPECFEM3D.f90.
% See also writeSTATIONS.m
%
% In Perl, this can be done as follows:
%   open(IN,${stations_seis}); @temp = <IN>; $nrec = $temp[0]; chomp($nrec);
%   print CSH "tail -n $nrec ${stations_seis} > stemp0 \n";
%   print CSH "awk '{print \$4,\$3}' stemp0 > temp \n";
%
% calls xxx
% called by xxx
%

function [rlon,rlat,relev,rburial,stnm,netwk] = readSTATIONS(filename)

% read in the STATIONS file
[stnm,netwk,rlat,rlon,relev,rburial] = textread(filename,'%s%s%f%f%f%f','headerlines',1);

%=======================================================================
