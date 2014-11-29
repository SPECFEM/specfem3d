%
% function 
% Carl Tape, 25-Jan-2010
%
% 
%
% calls wave2d_gll2cell.m, wave2d_cell2gll.m
% called by xxx
%

function m_out = wave2d_m_gll2cell(m_in,xg,yg,Xc,Yc,nxc,nyc,nmod_src,iopt)

nc = length(Xc(:));     % regular mesh of cells
ng = length(xg);        % irregular mesh GLL points

if iopt==1      % gll2cell
    nmod_out = nc+nmod_src;
    if length(m_in) ~= ng+nmod_src, error('check dimensions of input'); end
    m_out                = zeros(nmod_out,1);
    
    %iGLL                 = wave2d_gll2cell(xg,yg,xc,yc);
    %m_out(1:nc)          = m_in(iGLL);
    m_out(1:nc)          = wave2d_gll2cell(xg,yg,m_in(1:ng),Xc,Yc);
    
    m_out(nc+1:nmod_out) = m_in(ng+1:ng+nmod_src);
    
else            % cell2gll
    nmod_out = ng+nmod_src;
    if length(m_in) ~= nc+nmod_src, error('check dimensions of input'); end
    m_out                = zeros(nmod_out,1);
    
    %m_out(1:ng)          = wave2d_cell2gll(xg,yg,xc,yc,m_in(1:nc),nxc,nyc);
    M_in                 = reshape(m_in(1:nc),nyc,nxc);
    m_out(1:ng)          = wave2d_cell2gll(xg,yg,Xc,Yc,M_in);
    
    m_out(ng+1:nmod_out) = m_in(nc+1:nc+nmod_src);
end

%=========================================================