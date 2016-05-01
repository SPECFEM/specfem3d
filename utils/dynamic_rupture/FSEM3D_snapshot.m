%FSEM3D_SNAPSHOT reads fault data at a given time
%
% d = FSEM3D_snapshot(isnap, [dir, fault])
%
% INPUTS	isnap	snapshot index, as in Snapshot*.bin file names
%		dir	["."] directory containing the SPECFEM3D output data Snapshot*.bin
%		fault	[1] fault id
%
% OUTPUTS	d	structure containing fault fields:
%			X,Y,Z	  fault node coordinates (km)
%			Dx,Dz	  slip (m)
%			Vx,Vz	  slip rate (m/s)
%			Tx,Ty,Tz  stress change (MPa)
%			S	  slip "state" variable in friction law (m)
%			Sg	  strength relative to initial stress (MPa)
%			Trup	  rupture time (s)
%			Tpz	  process zone time, when slip=Dc (s)
%
% Jean-Paul Ampuero	ampuero@erdw.ethz.ch 
% 19/01/2011: modified by Percy Galvez percy.galvez@sed.ethz.ch 
% 02/03/2015: modified by Kangchen Bai kbai@caltech.edu
% WARNING : By default using double precision. Have to specify 'single' if
% otherwise.

function d = FSEM3D_snapshot(isnap,DATA_DIR,fault,precision)

NDAT = 14; 

if nargin<2, DATA_DIR = '.', precision = 'double'; end
if nargin<3, fault = 1, precision = 'double'; end
if nargin<4, precision = 'double'; end

len = 8;

if strcmp(precision,'single') , len = 4; end

BinFile = sprintf('%s/Snapshot%u_F%u.bin',DATA_DIR,isnap,fault);

if ~exist(BinFile,'file'), error(sprintf('File %s does not exist',BinFile)), end
fid=fopen(BinFile);
%BinRead = fread(fid,[1,inf],'single')' ;
%fclose(fid);
for ii = 1:1:NDAT
    number = fread(fid,1,'int');
    number
    N = number/len;
    BinRead(:,ii) = fread(fid,N,precision);
    number = fread(fid,1,'int');
end
fclose(fid);



d.X  = BinRead(:,1)/1e3; % in km
d.Y  = BinRead(:,2)/1e3; % in km
d.Z  = BinRead(:,3)/1e3; % in km
d.Dx = BinRead(:,4);
d.Dz = BinRead(:,5);
d.Vx = BinRead(:,6);
d.Vz = BinRead(:,7);
d.Tx = BinRead(:,8); % in MPa
d.Ty = BinRead(:,9);
d.Tz = BinRead(:,10); % in MPa
d.S  = BinRead(:,11);
d.Sg = BinRead(:,12); % in MPa
d.Trup = BinRead(:,13);
d.Tpz = BinRead(:,14);


return
