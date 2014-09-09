%
% wave2d_write_hessian.m
% Carl Tape, 01-Feb-2009
%
% 
% 
% calls xxx
% called by xxx
%

function wave2d_write_hessian(H,eids,odir,nprint)

disp(['wave2d_write_hessian.m: writing files to ' odir]);

[a,b] = size(H);
if a ~= b, error('H must b square'); end
nsrc = a;

nprintmax = (nsrc-1)*nsrc/2;
if nprint > nprintmax, error(sprintf('nprint (%i) must be less than nprintmax (%i)',nprint,nprintmax)); end

% construct the matrix of data-covariance normalization
i_ind = zeros(nsrc,nsrc);
j_ind = zeros(nsrc,nsrc);
for i = 1:nsrc
    for j = 1:nsrc
        i_ind(i,j) = i;
        j_ind(i,j) = j;
    end
end

% sort all entries of the Hessian -- note that some are negative
% we take the upper traingular elements, excluding the diagonal
Htri = triu(H,1); itri = triu(i_ind,1); jtri = triu(j_ind,1);
Hcols = [Htri(:) itri(:)  jtri(:) ];
bkeep = sum(Hcols,2)~=0; Hcols = Hcols(bkeep,:);
[Hsort,iHsort_good] = sortrows(Hcols, -1);
[Hsort,iHsort_bad]  = sortrows(Hcols, 1);

% make a list of the largest N positive and negative Hessian elements
filename = [odir 'hessian_elements_good'];
fid = fopen(filename,'w');
fprintf(fid,'Largest positive non-diagonal elements of Hessian:\n');
fprintf(fid,'%10s%6s%6s%12s%12s\n','Hij','i','j','eid-i','eid-j');
for ix = 1:nprint
    iH = iHsort_good(ix);
    fprintf(fid,'%10.2e%6i%6i%12s%12s\n',Hcols(iH,:),eids{Hcols(iH,2)},eids{Hcols(iH,3)});
end
fclose(fid);

filename = [odir 'hessian_elements_bad'];
fid = fopen(filename,'w');
fprintf(fid,'Largest negative non-diagonal elements of Hessian:\n');
fprintf(fid,'%10s%6s%6s%12s%12s\n','Hij','i','j','eid-i','eid-j');
for ix = 1:nprint
    iH = iHsort_bad(ix);
    fprintf(fid,'%10.2e%6i%6i%12s%12s\n',Hcols(iH,:),eids{Hcols(iH,2)},eids{Hcols(iH,3)});
end
fclose(fid);

%===================================================================
