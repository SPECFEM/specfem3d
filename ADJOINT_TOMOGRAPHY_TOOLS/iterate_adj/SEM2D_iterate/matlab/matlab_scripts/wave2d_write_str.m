%
% function wave2d_write_str(filename,x,y,kappa,mu,rho,B)
% Carl Tape, 25-Jan-2010
%
% This function writes out a structural model that can be read in by
% wave2d.f90.
%
% calls xxx
% called by xxx
%

function wave2d_write_str(filename,x,y,kappa,mu,rho,B)

disp(['writing to ' filename]);

% ! CURRENT MODEL (synthetics)
% write(19,'(6e20.10)') x_plot(iglob), z_plot(iglob), &
% kappa_syn(i,j,ispec), mu_syn(i,j,ispec), rho_syn(i,j,ispec), &
% log( beta_syn(i,j,ispec) / beta0 )
fid = fopen(filename,'w');
for ii = 1:length(x)
    fprintf(fid,'%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n',...
        x(ii), y(ii), kappa(ii), mu(ii), rho(ii), B(ii) );
end
fclose(fid);

%=========================================================