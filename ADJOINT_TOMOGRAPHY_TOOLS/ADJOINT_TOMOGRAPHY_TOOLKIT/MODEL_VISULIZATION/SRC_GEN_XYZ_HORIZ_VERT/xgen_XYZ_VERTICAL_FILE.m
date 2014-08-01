clear all
close all 
% cross section A, Carpathian-Adriatic slab, lon1=5(7.5),lat1=38(39); lon2=45(42.5),lat2=48(48); 
% cross section B, Hellenic arc, lon1=29, lat1=30; lon2=10, lat2=65;
% cross section C, Eifle Hotspot, lon1=-5, lat1=49;  lon2=30, lat2=47.5;
% cross section D, Calabria arc, lon1=-5, lat1=39;  lon2=30,lat2=37;
% cross section E, Alps,  lon1=10, lat1=30; lat2=10.1, lat2=65;
% cross section F, Apennies, lon1=-5(0),lat1=34(37); lon2=35(35), lat2=49(49); 
% cross section G, Scandiawia, lon1=10, lat1=68; lon2=15; lat2=48. 


% cross section A, Carpathian-Adriatic slab, lon1=15,lat1=42.1 ; lon2=35,lat2=47.1; 
% cross section B, Hellenic arc, lon1=29, lat1=32.5; lon2=19, lat2=52.5;
% cross section C, Alps,  lon1=10, lat1=35; lat2=10.1, lat2=50;
% cross section D, Scandiawia, lon1=11, lat1=65; lon2=14.5; lat2=50. 
% cross section E, Eifle Hotspot, lon1=0, lat1=50.5 ;  lon2=20, lat2=47.5 ;
% cross section F, Calabria arc, lon1=0, lat1=39.5;  lon2=20,lat2=38.5;
% cross section G, Apennies, lon1=0,lat1=32 ; lon2=20, lat2=47; 
% cross section H, Iaptus ocean, lon1=-10, lat1=61; lon2=10, lat2=48; 
% cross section SEIS, lon1=20.74,lat1=37.63; lon2=34.65, lat2=67.9; 
% cross section SINGLE, lon1=38.4273, lat1=38.3134; lon2=-6.6734, lat2=70.9573; 
% cross section CALA, lon1=0; lat1=38.8; lon2=20; lat2=38.8

% comparison between EU30 and Wortel & Spakman 
% cross section APEN, lon1=10; lat1=35.5; lon2=30; lat2=51 
% cross section APENNEW, lon1=0; lat1=38.2; lon2=20; lat2=41.2
% cross section SCAN, lon1=9.5; lat1=50; lon2=18; lat2=60
% cross section Turkey, 1, lon1=27, lat1=35; lon2=27, lat2=42
% cross section Turkey, 2, lon1=31, lat1=35; lon2=31, lat2=42
% cross section Turkey, 3, lon1=40, lat1=35; lon2=40, lat2=42
% cross section Turkey  4, lon1=26, lat1=39; lon2=44, lat2=39
% cross section Turkey  5, lon1=26, lat1=37; lon2=44, lat2=37
% cross section ICE1, lon1=-25; lat1=64.5; lon2=-12.5, lat2=64.5
% cross section ICE2, lon1=-21, lat1=63; lon2=-12.5, lat2=66.5
% cross section APLINE1, lon1=3.5,lat1=46; lon2=17.5,lat2=42.5
% cross section APLINE2, lon1=6.4,lat1=51.5;lon2=12.2,lat2=43.5
% cross section APLINE3, lon1=11,lat1=43.5;lon2=19,lat2=51

% cross section Iapetus, lon1=-3,lat1=59; lon2=8;lat2=48; 
% cross section Tornquist, lon1=6.5,lat1=49; lon2=20,lat2=61; 
% cross section England, lon1=1.5;lat1=48; lon2=-7,lat2=59;

lon1=1.5;
lat1=48;

lon2=-7;
lat2=59;

fnm='XYZ_FILE_VERT/VERT_SLICE_England.xyz';


fid=fopen(fnm,'w');
rtop=1;
rbot=0.80; % 1200km depth 

%npt1=800;
%npt2=500; 

npt1=300;
npt2=200;


dr=(rtop-rbot)/npt2; 

[ylat,xlon]=gcwaypts(lat2,lon2,lat1,lon1,npt1); %generate points along great circle 

for i = 1:npt1 

	phi=xlon(i)*pi/180.0;
	theta=(90.0-ylat(i))*pi/180.0;

	for j = 1:npt2
		r=rtop-(j-1)*dr;
		
		x=r*sin(theta)*cos(phi);
		y=r*sin(theta)*sin(phi);
		z=r*cos(theta);

		fprintf(fid,'%f  %f  %f\n',x,y,z);

	end 
end 
fclose(fid);
% cross section Alpine 1, lon1=3.5, lat1=46; lon2=16, lat2=43
% cross section Alpine 2, lon1=7.5, lat1=50.2; lon2=12.5, lat2=43;
% cross section Alpine 3, lon1=10.5, lat1=43; lon2=16.5, lat2=49


