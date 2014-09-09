%
% writeCMT_psmeca.m
% CARL TAPE, 02-Feb-2007
% printed xxx
%
% This inputs a set of moment tensors in units of N-m and outputs four
% different files for plotting in GMT using psmeca, which assumes units of
% dyne-cm.  The four files differ only in what label to use plotting above
% the beach balls.
%
% moment tensor M = [Mrr Mtt Mpp Mrt Mrp Mtp]
% From psmeca (-Sm) : mrr, mtt, mff, mrt, mrf, mtf in 10*exponent dynes-cm
%
% calls CMT2m0.m
% called by slab_readCMT.m, socal_quakes_2007.m
%

%function writeCMT_psmeca(filename,date,lat,lon,dep,M,eid,opts)
function writeCMT_psmeca(filename,date,lat,lon,dep,M,eid,isource)

ncmt = length(date);

% make sure M is N by 6
[a,b] = size(M); if b ~= 6, M = M'; end
[a,b] = size(M); if b ~= 6, error('dimension of M must be N by 6'); end

% if there is no isource input, then just default to all ones
if nargin < 8, isource = ones(ncmt,1); end

% convert moment tensor from N-m to dyne-cm
M = 1e7 * M;

% exponent for computing magnitude in psmeca
M0 = CMT2m0(1,M);
iexp_all = floor(log10(M0));

% for labeling the moment magnitude
Mw = m02mw(1,CMT2m0(1,M*1e-7));   % M0 must be in N-m

% options
%ilab = opts{1};     % type of label for beach ball

% write 6 different versions of the file, each having different labels
for ilab = 0:5
    if ilab==0, ext = ''; end
    if ilab==1, ext = '_eid'; end
    if ilab==2, ext = '_Mw'; end
    if ilab==3, ext = '_year'; end
    if ilab==4, ext = '_date'; end
    if ilab==5, ext = '_all'; end
    disp(['writeCMT_psmeca.m : extension is ' ext]);

    % write to file for GMT plotting
    file1 = [filename '_psmeca' ext];
    fid = fopen(file1,'w');
    for ii = 1:ncmt

        % title for beach ball
        switch ilab
            case 0, cmtlabel = '  ';
            case 1, cmtlabel = char(eid{ii});
            case 2, cmtlabel = sprintf('%.2f',Mw(ii));  
            case 3, cmtlabel = datestr(date(ii),10); 
            case 4, cmtlabel = datestr(date(ii),29); 
            case 5, cmtlabel = sprintf('%4i-M%.1f-%i',year(date(ii)),Mw(ii),isource(ii)); 
        end

        fac = 10^-iexp_all(ii);

        fprintf(fid,'%14.6f%14.6f%14.6f%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%4i%14.6e%14.6e%16s\n',...
            lon(ii), lat(ii), dep(ii),...
            M(ii,1)*fac, M(ii,2)*fac, M(ii,3)*fac,...
            M(ii,4)*fac, M(ii,5)*fac, M(ii,6)*fac,...
            iexp_all(ii),...
            lon(ii), lat(ii), cmtlabel);   
    end
    fclose(fid);
end

% write a list of event IDs
file2 = [filename '_eid'];
fid = fopen(file2,'w');
for ii=1:ncmt, fprintf(fid,'%s\n',char(eid{ii})); end
fclose(fid);

disp(' writing psmeca file for GMT plotting...');
disp([' output file : ' file1]);
disp([' output file : ' file2]);
disp([' number of CMT solutions : ' num2str(ncmt)]);

%======================================================
