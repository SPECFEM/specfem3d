%
% function 
% Carl Tape, 25-Jan-2010
%
% Given an input grid (xg,yg) for which a function is defined, this
% function returns the index of the nearest (xg,yg) point to each input
% (x,y) point.  This provides a quick vector for nearest neighbor
% interpolation, which is of course CRUDE.
%
% calls xxx
% called by xxx
%

%function [inds, dmins] = wave2d_gll2cell(xg,yg,xc,yc)
function Fc = wave2d_gll2cell(xg,yg,fg,Xc,Yc)

Fc = griddata(xg,yg,fg,Xc,Yc,'cubic');

% nc = length(xc);    % number of input points (regular mesh)
% ng = length(xg);    % number of GLL gridpoints (irregular mesh)
%
% % loop over input points
% inds = zeros(nc,1);
% dmin = zeros(nc,1);
% for ii=1:nc
%     dtemp = sqrt( (xc(ii)-xg).^2 + (yc(ii)-yg).^2 );
%     [dmin,imin] = min( dtemp );
%     inds(ii) = imin(1);
%     dmins(ii) = dmin(1);
% end

iplot = 0;
if iplot==1
    figure; plot(xc,yc,'k.',xg(inds),yg(inds),'ro');
    legend('cell center','closest GLL point');
    axis equal, xlabel('x'); ylabel('y');
    figure; plot(dmins,'k.'); axis tight;
    xlabel('Cell index'); ylabel('Distance between cell and GLL, meters');
end

%=========================================================