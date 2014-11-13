%
% function [xf,yf,mf,bf,rms] = linefit(x, y)
% CARL TAPE, 16-June-2005
% printed xxx
%
% This function fits data to a line.
% Note: With t points, rms = 0.
% EXAMPLE is shown below.
% 
% calls xxx
% called by c101D.m, c163E.m, c111C.m
% 

function [xf,yf,mf,bf,rms] = linefit(x,y)

% column vectors
x = x(:); y = y(:);

% remove NaN or inf
igood = find( ~isnan(x) .* ~isinf(x) .* ~isnan(y) .* ~isinf(y) == 1 );
x = x(igood);
y = y(igood);

% sample xpts
xf = linspace(min(x), max(x), 100)';

% least-squares fitting
P = polyfit(x,y,1);
yf = polyval(P,xf);

mf = P(1);
bf = P(2);
res = polyval(P,x) - y;
rms = sqrt( (res' * res) / length(res) );

if 0==1
    % EXAMPLE (c101D.m)
    Srvec = [28.8 48.8 62.3 66.6 72.5 78.6];
    SrRatio = [0.7151 0.7130 0.7123 0.7112 0.7113 0.7111];
    [xf,yf,mf,bf,rms] = linefit(SrRatio, 1./Srvec);
    figure; plot(xf,yf,'r--',SrRatio,1./Srvec,'b.');
    xlabel('^{87}Sr / ^{86}Sr');
    ylabel('1/[Sr] (1/ppb)');
end

%-----------------------------------------------------------
% older version which requires OPTIMIZATION toolbox

% % initial guess based on the first and last points in data
% mg = (y(end) - y(1)) / (x(end) - x(1));
% bg = y(1) - mg*x(1);
% 
% ydata = mg*x + bg;
% fun = inline('x(1)*xdata + x(2)', 'x', 'xdata');
% F = lsqcurvefit(fun, [x(1) x(end)], x, y);
%
% mf = F(1);
% bf = F(2);
% yf = mf*xf + bf;

%===============================================================
