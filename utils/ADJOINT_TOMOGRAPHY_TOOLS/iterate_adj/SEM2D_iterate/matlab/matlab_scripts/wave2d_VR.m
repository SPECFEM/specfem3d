%
% function wave2d_VR
% Carl Tape, 25-Jan-2010
%
% 
%
% calls xxx
% called by wave2d_cg_run.m
%

function [VR_i, VR_s, VR_tot] = wave2d_VR(DT1,DT2,sigmaDT,Cd,Dinds)

% D_inds is nsrc x 2
nsrc = length(Dinds);
n = length(DT1);

chi1 = DT1.^2 ./ Cd;
chi2 = DT2.^2 ./ Cd;

VR_i = log( chi1 ./ chi2 );

% find indices of measurements for which both DT1 and DT2 are less than sigmaDT
% --> VR = 0 for these
ivr = find( and( abs(DT1) < sigmaDT, abs(DT2) < sigmaDT ) );
disp(sprintf('%i / %i measurements with DT < %.2f',length(ivr),n,sigmaDT));
%for ii=1:length(ivr)
%    disp(sprintf('%10i%10.4f%10.4f',ivr(ii),DT1(ivr(ii)),DT2(ivr(ii))))
%end
VR_i(ivr) = 0;

VR_s = zeros(nsrc,1);
for ii=1:nsrc
    inds = [Dinds(ii,1) : Dinds(ii,2)];
    VR_s(ii) = sum( VR_i(inds) );
end

VR_tot = log( sum(chi1) / sum(chi2) );

if 1==1
    figure; nr=3; nc=1;

    subplot(nr,nc,1); hold on;
    plot( DT1, DT2, '.');
    xlabel('DT - model 1'); ylabel('DT - model 2');
    axis equal, axis tight, grid on;
    ax1 = axis; val = min(ax1([2 4]));
    plot([-1 1]*val,[-1 1]*val,'r--','linewidth',3);

    subplot(nr,nc,2); hold on;
    plot( DT2 - DT1, '.'); grid on; xlim([0 n])
    xlabel('index'); ylabel('DT2 - DT1')

    subplot(nr,nc,3); hold on;
    plot( VR_i, '.'); grid on; xlim([0 n])
    xlabel('index'); ylabel('VR')
    
    orient tall, wysiwyg
end

%=========================================================
