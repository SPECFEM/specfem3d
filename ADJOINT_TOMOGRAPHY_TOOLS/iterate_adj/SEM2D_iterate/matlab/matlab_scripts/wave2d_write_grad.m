%
% function 
% Carl Tape, 25-Jan-2010
%
% 
%
% calls xxx
% called by xxx
%

function wave2d_write_grad(filename,gk,pk)

fid = fopen(filename,'w');
for ii = 1:length(gk)
    fprintf(fid,'%20.10e%20.10e\n',gk(ii), pk(ii) );
end
fclose(fid);

%=========================================================