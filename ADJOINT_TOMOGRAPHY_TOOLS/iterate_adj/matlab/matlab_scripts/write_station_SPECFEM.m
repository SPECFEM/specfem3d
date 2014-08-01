%
% function write_station_SPECFEM(filename,rlon,rlat,relev,rburial,stnm,netwk)
% CARL TAPE, 15-March-2007
% printed xxx
%
% This function inputs a set of values and names for stations, and it
% outputs a stations file suitable for:
%    (1) plotting in GMT
%    (2) running SPECFEM3D
%
% In Perl, this can be done as follows:
%   open(IN,${stations_seis}); @temp = <IN>; $nrec = $temp[0]; chomp($nrec);
%   print CSH "tail -n $nrec ${stations_seis} > stemp0 \n";
%   print CSH "awk '{print \$4,\$3}' stemp0 > stemp \n";
%   print CSH "psxy stemp ${station_info2} $J $R -K -V -O >> $psfile \n";
%
% calls xxx
% called by xxx
%

function write_station_SPECFEM(filename,rlon,rlat,relev,rburial,stnm,netwk)

% number of receivers
nrec = length(rlat);

% GMT format
ofile = [filename '_gmt'];
disp([' Writing the stations file ' ofile ]);
fid = fopen(ofile,'w');
for ii=1:nrec
    fprintf(fid,'%14.6f%14.6f%7s%4s%12.2f%12.2f\n',...
        rlon(ii),rlat(ii),stnm{ii},netwk{ii},relev(ii),rburial(ii));   
end
fclose(fid);

% SPECFEM format
ofile = [filename '_specfem'];
disp([' Writing the stations file ' ofile ]);
fid = fopen(ofile,'w');
fprintf(fid,'%10i\n',nrec);
for ii=1:nrec
    fprintf(fid,'%7s%4s%14.6f%14.6f%12.2f%12.2f\n',...
        stnm{ii},netwk{ii},rlat(ii),rlon(ii),relev(ii),rburial(ii));   
end
fclose(fid);
    
%=======================================================================
