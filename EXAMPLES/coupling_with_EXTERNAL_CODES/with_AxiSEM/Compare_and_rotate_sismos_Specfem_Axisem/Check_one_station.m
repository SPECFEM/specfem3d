SIS='./Sismo_axisem'; % axisem sismogtram directory 
prefix='S0001_SY_disp_post_mij_conv0000_'; % axisem prefix name file
Specfem_rot=load('Sismo_specfem/SY_0001_ENZ.displ');     % specfem rotated sesmograms 

T0=710; % origin time chooen in reformat.par


%


CN=[prefix,'N.dat'];
CE=[prefix,'E.dat'];
CZ=[prefix,'Z.dat'];


fichier=fullfile(SIS,CN);
displ_N_axiSEM=load(fichier);

fichier=fullfile(SIS,CE);
displ_E_axiSEM=load(fichier);

fichier=fullfile(SIS,CZ);
displ_Z_axiSEM=load(fichier);

%
figure; 

%E
subplot(3,1,1)
plot(Specfem_rot(:,1)+T0-Specfem_rot(1,1),Specfem_rot(:,2))
hold on;
plot(displ_E_axiSEM(:,1)-displ_E_axiSEM(1,1),displ_E_axiSEM(:,2),'r')
legend('Hybrid', 'AxiSEM');
ylabel('E')
%N
subplot(3,1,2)
plot(Specfem_rot(:,1)+T0-Specfem_rot(1,1),Specfem_rot(:,3))
hold on;
plot(displ_N_axiSEM(:,1)-displ_N_axiSEM(1,1),displ_N_axiSEM(:,2),'r')
ylabel('N')

%Z
subplot(3,1,3)
plot(Specfem_rot(:,1)+T0-Specfem_rot(1,1),Specfem_rot(:,4))
hold on;
plot(displ_Z_axiSEM(:,1)-displ_Z_axiSEM(1,1),displ_Z_axiSEM(:,2),'r')
xlabel('time (s)')
ylabel('Z')

linkaxes

xlim([710 1800])