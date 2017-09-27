%FSEM3D_moment_area computes seismic moment or potency and rupture area
%
% [M,A] = FSEM3D_moment_area(ngll,isnap, [data_dir,db_dir,Dmin])
%
% INPUTS
% ngll 	number of GLL points per element edge (parameter NGLL of SPECFEM3D)
% isnap	last snapshot index in Snapshot*.bin file names (contains final slip)
% dat_dir	["OUTPUT_FILES"] directory containing the files Snapshot*.bin
% db_dir	["OUTPUT_FILES/DATABASES_MPI"] directory containing the files proc*fault_db.bin
% Dmin	  [1e-3] minimum slip to define the rupture area
% s_or_d  ['single'] Precision of SPECFEM3D outputs, 'single' or 'double'
% mu_type  [1]  Type of shear modulus setting (see nested function my_mu):
%     1 = sets mu=1 and outputs Px and Pz are potency
%     2 = an ad hoc user-defined setting (use it as an example to implement your own)
%
% OUTPUTS	
% Px	seismic potency or moment along-strike (integral of along-strike slip)
% Pz	seismic potency or moment along-dip (integral of along-dip slip)
% A	area of regions with slip > Dmin
%
% NOTES 
% If multiple faults exist, this function only processes the first fault
%   
% Jean-Paul Ampuero	ampuero@gps.caltech.edu

%		fault	[1] fault id

function [Px,Pz,A] = FSEM3D_moment_area(ngll,isnap,data_dir,db_dir,Dmin,s_or_d,mu_type) %,fault)

if nargin<3, data_dir = 'OUTPUT_FILES'; end
if nargin<4, db_dir = 'OUTPUT_FILES/DATABASES_MPI'; end
if nargin<5, Dmin=1e-3; end 
if nargin<6, s_or_d='single'; end 
if nargin<7, mu_type=1; end

%if nargin<6, fault = 1; end
fault=1;

% read slip
dat = FSEM3D_snapshot(isnap,data_dir,fault,s_or_d);
Dx=dat.Dx;
Dz=dat.Dz;
Z = dat.Z;
clear dat

list = dir([db_dir '/*fault_db.bin']); 
nproc = length(list);

Px=0;
Pz=0;
A=0;
i0=0;

for p=1:nproc

  fid=fopen([db_dir '/' list(p).name]);
  dat = fread(fid,3,'int');
 % WARNING: reading the first fault only
  dat = fread(fid,4,'int');
  nspec=dat(2);
  nglob=dat(3);
  if (nspec==0), fclose(fid); continue; end % if this proc does not contain the fault, skip to next proc
  dat = fread(fid,nspec*ngll^2+2,'int');
  ibool = dat(2:end-1);
  dat = fread(fid,1,'int');
  jacw = fread(fid,nspec*ngll^2,s_or_d)';
  fclose(fid);

  Dpx = Dx(i0+ibool);
  Dpz = Dz(i0+ibool);
 
 % compute area
  Dp2 = Dpx.*Dpx+Dpz.*Dpz;
  rup_mask = (Dp2>Dmin^2);
  A = A + jacw*rup_mask;
 
 % compute potency or moment
  mu = my_mu(mu_type);
  mu = mu(:).*rup_mask(:);
  Dpx = Dpx.*mu;
  Dpz = Dpz.*mu;
  
  Px = Px + jacw*Dpx;
  Pz = Pz + jacw*Dpz;
  
  i0 = i0+nglob;
end

%---
% This is a nested function: it can use and modify all variables of the function above
function mu=my_mu(mu_type)

switch mu_type

  case 1
    mu = 1;
  
  case 2
    mu=3200*3200*2110;
    if (Z(i0+ibool)<-4) 
      mu = 3200*3200*2720;
    end
    if (Z(i0+ibool)<-24) 
      mu = 3700*3700*2790; 
    end
    if(Z(i0+ibool)<-46)
      mu = 4550*4550*3380;
    end

  end

end

%---
end
