%
% function wave2d_tsvd_hessian(H,dnorm)
% Carl Tape, 25-Jan-2010
%
% This function examines possible regularization of H via truncated
% singular value decomposition.
%
% calls tsvd.m
% called by wave2d_cg_run.m
%

function [mu_all] = wave2d_tsvd_hessian(H,dnorm)

nsrc = length(dnorm);

% Matlab SVD
[~,S,~] = svd(H);
s = diag(S);

pinds = [1:nsrc]';
[mu_all,rss,f_r_ss] = tsvd(dnorm,H,pinds);   % KEY: TSVD

% norms of mu vectors
mu_norm = zeros(nsrc,1);
for ip = 1:nsrc
    mu_norm(ip) = norm(mu_all(:,ip));
end

figure; nr=2; nc=2;
xlab1 = 'p, singular value index';
xlab2 = 'p, singular value truncation index';
ylab1 = 'singular value';
ylab2 = 'misfit : dot[ d - H*mu(p), d - H*mu(p) ]';

subplot(nr,nc,1); plot(pinds,s,'.','markersize',20);
grid on; xlabel(xlab1); ylabel(ylab1);
subplot(nr,nc,2); semilogy(pinds,s,'.','markersize',20);
grid on; xlabel(xlab1); ylabel(ylab1);
subplot(nr,nc,3); plot(pinds,rss,'.','markersize',20);
grid on; xlabel(xlab2); ylabel(ylab2);
subplot(nr,nc,4); semilogy(pinds,rss,'.','markersize',20);
grid on; xlabel(xlab2); ylabel(ylab2);
orient tall, wysiwyg

figure; nr=2; nc=2;
ylab3 = 'norm of mu vector';

subplot(nr,nc,1); semilogy(pinds,s,'.','markersize',20);
grid on; xlabel(xlab1); ylabel(ylab1);
subplot(nr,nc,2); semilogy(pinds,rss,'.','markersize',20);
grid on; xlabel(xlab2); ylabel(ylab2);
subplot(nr,nc,3); semilogy(pinds,mu_norm,'.','markersize',20);
grid on; xlabel(xlab2); ylabel(ylab3);
subplot(nr,nc,4); loglog(mu_norm,rss,'.-','markersize',20);
grid on; xlabel(ylab3); ylabel(ylab2);
orient tall, wysiwyg

% mu vectors
figure; nr=2; nc=1;
subplot(nr,nc,1); plot(pinds,mu_all,'.-');
grid on; xlabel('event index'); ylabel('elements of mu vectors');
xlim([0 nsrc+1]);
subplot(nr,nc,2); plot(mu_all','.-');
grid on; xlabel(xlab2); ylabel('elements of mu vectors');
xlim([0 nsrc+1]);
orient tall, wysiwyg

%=========================================================
