%FSEM3D_moment_area computes seismic moment and rupture area
%
% [M,A] = FSEM3D_moment_area(ngll,isnap, [data_dir,db_dir,Dmin])
%
% INPUTS	ngll 	number of GLL points per element edge (parameter NGLL of SPECFEM3D)
%		isnap	last snapshot index in Snapshot*.bin file names (contains final slip)
%		dat_dir	["OUTPUT_FILES"] directory containing the files Snapshot*.bin
%		db_dir	["OUTPUT_FILES/DATABASES_MPI"] directory containing the files proc*fault_db.bin
%		Dmin	minimum slip to define the rupture area
%
% OUTPUTS	Px	seismic potency along-strike (integral of along-strike slip)
% 		Pz	seismic potency along-dip (integral of along-dip slip)
%		A	area of regions with slip > Dmin
%
% WARNING: implemented only to process the first fault
% 
% Jean-Paul Ampuero	ampuero@gps.caltech.edu

%		fault	[1] fault id

function [Px,Pz,A] = FSEM3D_moment_area(ngll,isnap,data_dir,db_dir,Dmin) %,fault)

if nargin<3, data_dir = 'OUTPUT_FILES'; end
if nargin<4, db_dir = 'OUTPUT_FILES/DATABASES_MPI'; end
if nargin<5, Dmin=1e-3; end 

%if nargin<6, fault = 1; end
fault=1;

% read slip
dat = FSEM3D_snapshot(isnap,data_dir,fault);
Dx=dat.Dx;
Dz=dat.Dz;
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
  dat = fread(fid,nspec*ngll^2+2,'float');
  jacw = dat(2:end-1)';
  fclose(fid);

  Dpx = Dx(i0+ibool);
  Dpz = Dz(i0+ibool);
  Px = Px + jacw*Dpx;
  Pz = Pz + jacw*Dpz;
  Dp2 = Dpx.*Dpx+Dpz.*Dpz;
  A = A + jacw*(Dp2>Dmin^2);

  i0 = i0+nglob;

end
