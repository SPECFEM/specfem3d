
clear all;
close all;

%Output file located in ~/SPECFEM3D/OUTPUT_FILES
BinFile = 'moviedata033000'; 
par = 3;
scale = 3;
SLIM = [-scale scale];

NDAT = 6 ;
fid = fopen(BinFile);
BinRead = fread(fid,[1,inf],'single')';
fclose(fid);

BinRead = reshape( BinRead(:),length(BinRead)/(NDAT),NDAT);
BinRead = BinRead(2:end-1,:);

        d.X  = BinRead(:,1)/1e3; % in km
        d.Y  = BinRead(:,2)/1e3; % in km
        d.Z  = BinRead(:,3)/1e3; % in km
        d.dX = BinRead(:,4);
        d.dY = BinRead(:,5);
        d.dZ = BinRead(:,6);

switch 1
    case(par == 1)
        d.dD  = d.dX ;
    case(par == 2)
        d.dD  = d.dY ;
    case(par == 3)
        d.dD  = d.dZ ;
end  
        
XLIM = [0 120];
YLIM = [0 300];
% Shift to Wendt's coordinate system
d.Y = 180-d.Y;
d.X = d.X;

ind = find(d.Z == 0);

dxsurf  = d.Y(ind);
dysurf  = d.X(ind);
ddDsurf = d.dD(ind);

indrange =find((dxsurf > 0 & dxsurf < 120) & (dysurf > -50 & dysurf < 300));

x  = dxsurf(indrange);
y  = dysurf(indrange);
dz = ddDsurf(indrange);

scatter3(x,y,dz,[],dz,'.');
axis('tight')
xlabel('Along strike (km)');
ylabel('Down dip (km)');

zlim([-1 1]);
caxis([-scale scale]);

az= 55 ; el=70;
view([az el]);

colorbar

