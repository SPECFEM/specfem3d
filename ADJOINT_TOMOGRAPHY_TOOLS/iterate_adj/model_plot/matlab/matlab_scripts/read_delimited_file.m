%
% [data,ilabs,labels] = read_delimited_file(filename,ncol);
% CARL TAPE, 22-August-2007
% printed xxx
%
% This function reads in a delimited file.  It was designed to read files
% that were output from GMT using the pscoast command.
%
% NOTE: THERE MUST BE A DELIMITER TO START AND END THE FILE.
%
% calls xxx
% called by xxx
%

function [data,ilabs,labels,inds] = read_delimited_file(filename,ncol)

disp([' reading ' filename '...']);

% read all lines
lines = textread(filename,'%s','delimiter','\n','whitespace','');
nlines = length(lines);

% read in the LUR file and string delimiters
jj = 0; kk = 0;
data = NaN * ones(nlines,ncol);
for ii = 1:nlines
    tline = lines(ii);
    temp = str2num(char(tline));

    if ~isempty(temp)
        jj = jj + 1;
        data(ii,:) = temp;
    else
        kk = kk + 1;
        labels(kk) = tline;
        ilabs(kk) = ii;
    end
end
ilabs = ilabs(:);
labels = labels(:);

% indices for each segment
inds = [ilabs(1:end-1)+1 ilabs(2:end)-1];

%=========================================================================
