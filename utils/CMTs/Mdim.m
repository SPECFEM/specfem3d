function [Mout,n] = Mdim(M)
%MDIM checks to make sure that the matrix M is 6 x n; if not, the transpose is applied.

[a,b] = size(M);

if and(a~=6, b~=6)
    M
    disp(sprintf('M is %i by %i',a,b));
    error('neither dimension of M is 6');
end

if and(a==6, b==6)
    disp('Mdim.m: WARNING input M is 6 x 6');
    disp('        make sure that elements are properly ordered');
    
else
    % take transpose
    if a ~= 6
        M = M';
        disp(sprintf('Mdim.m: M is %i x %i',a,b));
        disp('Mdim.m: WARNING applying transpose to M');
    end
end

Mout = M;
[~,n] = size(Mout);
