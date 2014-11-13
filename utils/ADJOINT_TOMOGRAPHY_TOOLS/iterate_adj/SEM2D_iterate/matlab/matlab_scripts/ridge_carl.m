%
% function 
% Carl Tape (Tapio Schneider, ACM 118)
% 06-Nov-2006
% 
% This function inputs a design matrix, a data vector, and a vector of
% regularization parameters, and it returns three different curves that may
% be used to select the best parameter:
%    (1) L-curve and curvature
%    (2) generalized cross-validation curve (GCV)
%    (3) ordinary cross-validation (OCV), also known as 'leave-one-out' CV
%
% It is best to input a large number of regularization parameters, so that
% the min and max of the respective functions can be easily obtained.
%
% This program is copied in part from ridge_tapio.m
%
% NOTE THE PLOTTING OPTIONS.
% 
%-------------------------------------------------
% RIDGE  Ridge regression estimates.
%
%    Given a vector g, a design matrix X, and
%    a regularization parameter h, 
%           
%           [m, rss, mss, dof] = ridge_tapio(g, X, h) 
% 
%    returns the ridge regression estimate of the vector f in the
%    linear regression model
% 
%           g = X*f + noise.
%
%    Also returned are the residual sum of squares rss, the sum of
%    squares mss of the elements of the ridge regression estimate
%    m (the squared norm of m), and the effective number of
%    residual degrees of freedom dof.
%
%    If h is a vector of regularization parameters, the i-th column
%    m(:,i) is the ridge regression estimate for the regularization
%    parameter h(i); the i-th elements of rss and mss are the
%    associated residual sum of squares and estimate sum of squares.
%
%    If no regularization parameter h is given, generalized
%    cross-validation is used to determine the regularization
%    parameter. The chosen regularization parameter h and the value of
%    the GCV function are then returned as the fifth and sixth
%    output arguments 
%
%            [m, rss, mss, dof, h, G] = ridge_tapio(g, X);
% 
%    Adapted from various routines in Per Christian Hansen's
%    Regularization Toolbox.
%
% calls gcvfctn.m, curvature.m
% called by xxx
% 

function [m, rss, mss, Gvec, Fvec, dof, kap, iL, iGCV, iOCV] = ridge_carl(dvec, X, hvec)

% Size of inputs
[n, p]      = size(X); 
q           = min(n, p);
nh          = length(hvec);
if (min(hvec) < 0)
    error('Impossible regularization parameter h.')
end

% Initialize outputs
m         = zeros(p, nh);
rss         = zeros(nh, 1); 
mss      = zeros(nh, 1);
dof         = zeros(nh, 1);

% Compute SVD of X
[U, S, V]   = svd(X, 0);  
s           = diag(S);      % vector of singular values
s2          = s.^2;

% Coefficients in expansion of solution in terms of right singular vectors
fc          = U(:, 1:q)'*dvec;
zeta        = s .* fc;

% Treat each regularization parameter separately
for j = 1:nh    
    m(:, j) = V(:, 1:q) * (zeta ./ (s2 + hvec(j)^2));
    mss(j)  = sum(m(:, j).^2);
    rss(j)  = hvec(j)^4 * sum(fc.^2 ./ (s2 + hvec(j)^2).^2);
    dof(j)  = n - sum(s2./(s2 + hvec(j)^2));
end

% In overdetermined case, add rss of least-squares problem
if (n > p)
    rss = rss + sum((dvec - U(:, 1:q)*fc).^2);
end

%-----------------------
% determine the Lcurve pick (max curvature)

x1 = log10(rss);
y1 = log10(mss);

% % smooth curvature interpolation to get h_L
% num = 1000;
% xsmooth = linspace(x1(1),x1(end),1000);
% ysmooth = interp1(x1,y1,xsmooth,'cubic');
% [i0,kap_smooth] = curvature(xsmooth,ysmooth);
% rss_L = 10^xsmooth(i0);
% mss_L = 10^ysmooth(i0);
% h_L = 10^interp1(x1,log10(hvec),xsmooth(i0),'cubic');

% curvature, based on input h values alone
[iL,kap] = curvature(x1,y1);
%h_L = hvec(iL);
%rss_L = 10^x1(iL);
%mss_L = 10^y1(iL);

%-----------------------
% obtain GCV `best' solution and GCV curve

% GCV minimum -- 'exact' in the sense a minimization method is used
% [hmin, Gmin] = gcv_tapio(U, s, dvec, 'ridge');

% GCV minimum -- 'crude' in the sense that we coarsely sample the function
dof0 = n-q;
rss0 = sum((dvec - U(:, 1:q)*fc).^2);
Gvec = zeros(nh,1);
for j = 1:nh 
    Gvec(j) = gcvfctn(hvec(j), s2, fc, rss0, dof0);
end
[Gmin,iGCV] = min(Gvec);
%hmin = hvec(iGCV);

% GCV best model and L-curve point for h_GCV (= hmin)
%mod_min = inv(X'*X + hmin^2*eye(p))*X'*dvec;
%res = X*mod_min - dvec;
%rss_min = sum(res.^2);
%mss_min = sum(mod_min.^2);

% compute G for the Lcurve pick
%G_L = gcvfctn(h_L, s2, fc, rss0, dof0);

%-----------------------
% ordinary (leave-one-out) cross-validation

Fvec = ocv_carl(dvec, X, hvec);
[Fmin,iOCV] = min(Fvec);

%======================================================
% PLOTTING

lamL = hvec(iL);   GL = Gvec(iL);   rssL = rss(iL);   mssL = mss(iL);   kapL = kap(iL);    FL = Fvec(iL); 
lamF = hvec(iOCV); GF = Gvec(iOCV); rssF = rss(iOCV); mssF = mss(iOCV); kapF = kap(iOCV);  FF = Fvec(iOCV);
lamG = hvec(iGCV); GG = Gvec(iGCV); rssG = rss(iGCV); mssG = mss(iGCV); kapG = kap(iGCV);  FG = Fvec(iGCV); 

x1 = log10(rss);
y1 = log10(mss);
x2 = log10(hvec);
y2 = kap;
x3 = log10(hvec);
y3 = log10(Fvec);
x4 = log10(hvec);
y4 = log10(Gvec);

stx1 = ' Misfit norm, log10 RSS';
sty1 = ' Model norm, log10 MSS';
stx2 = ' Regularization parameter, log10 \lambda';
sty2 = ' Curvature of L-curve, \kappa(\lambda)';
stx3 = stx2;
sty3 = ' OCV function, log10 F(\lambda)';
stx4 = stx2;
sty4 = ' GCV function, log10 G(\lambda)';

%stfm = '%.4f';
stfm = '%.2e';
stlam_L   = [' \lambda-L = ' num2str(sprintf(stfm, lamL))];
stlam_ocv = [' \lambda-OCV = ' num2str(sprintf(stfm, lamF))];
stlam_gcv = [' \lambda-GCV = ' num2str(sprintf(stfm, lamG))];

%------------------------
figure; nr=2; nc=2;
msize = 8;
nlab = 10; ilabs = round(linspace(1,nh,nlab));
    
subplot(nr,nc,1); hold on;
plot(x1,y1,'.');
plot(log10(rssL),log10(mssL),'ko','markersize',8,'MarkerFaceColor','r');
plot(log10(rssG),log10(mssG),'kV','markersize',8,'MarkerFaceColor','g');
plot(log10(rssF),log10(mssF),'k^','markersize',8,'MarkerFaceColor','c');
axis tight; ax1 = axis; axis(axes_expand(ax1,1.1,1));
axy = axis; dx = axy(2)-axy(1);
for kk=1:nlab
   ii = ilabs(kk);
   text(x1(ii)+dx*0.02,y1(ii),[num2str(sprintf(stfm, hvec(ii)))],'fontsize',8,'color','b');
end
legend(' \lambda',stlam_L,stlam_gcv,stlam_ocv);
xlabel(stx1); ylabel(sty1); grid on;

subplot(nr,nc,2); hold on;
plot(x2,y2,'.');
plot(log10(lamL),kapL,'ko','markersize',8,'MarkerFaceColor','r');
plot(log10(lamG),kapG,'kV','markersize',8,'MarkerFaceColor','g');
plot(log10(lamF),kapF,'k^','markersize',8,'MarkerFaceColor','c');
axis tight; ax1 = axis; axis(axes_expand(ax1,1.1,1));
legend(' \lambda',stlam_L,stlam_gcv,stlam_ocv,'location','northwest');
xlabel(stx2); ylabel(sty2); grid on;

subplot(nr,nc,3); hold on;
plot(x3,y3,'.');
plot(log10(lamL),log10(FL),'ko','markersize',8,'MarkerFaceColor','r');
plot(log10(lamG),log10(FG),'kV','markersize',8,'MarkerFaceColor','g');
plot(log10(lamF),log10(FF),'k^','markersize',8,'MarkerFaceColor','c');
axis tight; ax1 = axis; axis(axes_expand(ax1,1.1,1));
legend(' \lambda',stlam_L,stlam_gcv,stlam_ocv,'location','northeast');
xlabel(stx3); ylabel(sty3); grid on;

subplot(nr,nc,4); hold on;
plot(x4,y4,'.');
plot(log10(lamL),log10(GL),'ko','markersize',8,'MarkerFaceColor','r');
plot(log10(lamG),log10(GG),'kV','markersize',8,'MarkerFaceColor','g');
plot(log10(lamF),log10(GF),'k^','markersize',8,'MarkerFaceColor','c');
axis tight; ax1 = axis; axis(axes_expand(ax1,1.1,1));
legend(' \lambda',stlam_L,stlam_gcv,stlam_ocv,'location','northwest');
xlabel(stx4); ylabel(sty4); grid on;

orient tall, wysiwyg, fontsize(9)
   
%======================================================
