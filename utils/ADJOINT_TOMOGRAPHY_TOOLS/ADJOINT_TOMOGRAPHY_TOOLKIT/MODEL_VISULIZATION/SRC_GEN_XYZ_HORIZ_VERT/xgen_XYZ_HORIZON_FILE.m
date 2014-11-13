clear all
close all 

depth=400;    % input depth km 

rnew=1.0-depth/6371.0;

fnm_output=['XYZ_FILE_HORIZ/DEPTH_SLICE_',num2str(depth,'%3.3i'),'.xyz'];


nglob=921600;  %total number of GLL points
fnm_input='SURFACE_GLL_POINTS.xyz';

fid_input=fopen(fnm_input,'r');
fid_output=fopen(fnm_output,'w');


for i = 1:nglob
	i
	array=fscanf(fid_input,'%f',2);
	lon=array(1);
	lat=array(2);


	phi=lon*pi/180.0;
	theta=(90.0-lat)*pi/180.0;

	xnew=rnew*sin(theta)*cos(phi);
	ynew=rnew*sin(theta)*sin(phi);
	znew=rnew*cos(theta);
	
	fprintf(fid_output,'%f		%f	  %f\n',xnew,ynew,znew);

end 

fclose(fid_input);
fclose(fid_output);
