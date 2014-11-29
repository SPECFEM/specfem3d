%
% function diff_files(file1,file2)
% Carl Tape, 25-Jan-2010
%
% This function reads in two files of numeric data, computes the norms of each column,
% then compares.
%
% calls xxx
% called by xxx
%

function diff_files(file1,file2)

d1 = load(file1); [nr1,nc1] = size(d1);
d2 = load(file2); [nr2,nc2] = size(d2);

if nr1 ~= nr2
    disp(sprintf('file 1 has %i rows',nr1));
    disp(sprintf('file 2 has %i rows',nr2));
    error('number of rows does not match');
end
if nc1 ~= nc2
    disp(sprintf('file 1 has %i columns',nc1));
    disp(sprintf('file 2 has %i columns',nc2));
    error('number of columns does not match');
end

format long

for ii=1:nc1
    disp(sprintf('column %2i%16.8e%16.8e%16.8e',...
        ii,norm(d1(:,ii)), norm(d2(:,ii)), norm( d1(:,ii) - d2(:,ii) ) / norm(d1(:,ii)) ))
end

if 0==1
    bdir = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_OUTPUT/';
    file1 = [bdir 'run_1500/READ_IN_CGF90_CDIAG/model_m0001/src_syn_m0001.dat'];
    file2 = [bdir 'run_1500/READ_IN/model_m0001/src_syn_m0001.dat'];
    diff_files(file1,file2)
end

%=========================================================