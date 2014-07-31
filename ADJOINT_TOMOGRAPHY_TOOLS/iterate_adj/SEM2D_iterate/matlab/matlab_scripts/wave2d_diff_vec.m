%
% function 
% Carl Tape, 25-Jan-2010
%
% 
%
% calls xxx
% called by xxx
%

function wave2d_diff_vec(m1,m2,m_inds,m_labs,ifig)

disp('wave2d_diff_vec: comparing two vectors:');

for ii=1:4
    inds = m_inds(ii,1):m_inds(ii,2);
    disp(sprintf('%20s%10.3e%10.3e%10.3e%10.3e',m_labs{ii},...
        norm(m1(inds)), norm(m2(inds)),...
        norm(m1(inds)-m2(inds)),...
        norm(m1(inds)-m2(inds)) / norm(m1(inds)) ))
    
    if ifig==1
        figure; nr=3; nc=1;
        subplot(nr,nc,1); plot( m1(inds), '.'); title(stpars{ii});
        subplot(nr,nc,2); plot( m2(inds), '.')
        subplot(nr,nc,3); plot( m1(inds), m2(inds), '.')
    end
end


%=========================================================