function plot_adj(data,syn,new1,new2,mtm_adj,ctp_adj,b,e)

subplot(2,1,1)
title('data, syn, recon_syn')
plot(data(:,1),data(:,2),'-k')
hold on
plot(syn(:,1),syn(:,2),'--r')
plot(new1(:,1),new1(:,2),'-','Color',[0.1210 0.6050 0])
plot(new2(:,1),new2(:,2),'-b');
xlim([b,e]);
xlabel('Time [s]')
legend('Data','Synthetic','MTM Recon_Syn', 'STP Recon_Syn')

subplot(2,1,2)
plot(mtm_adj(:,1),mtm_adj(:,2),'-','Color',[0.1210 0.6050 0]);
hold on;
plot(ctp_adj(:,1),ctp_adj(:,2),'-b');
xlim([b,e]);
legend('MT adjoint source', 'CC adjoint source');
