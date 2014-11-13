%
% socal_tomo_plots.m
% Carl Tape, 27-March-2009
%
% This file generates xyz points for which we want to extract an
% interpolation value from a SPECFEM volumetric field.  Once we have the
% function value, we can plot the cross-section in GMT.
%
% calls merge_arrays.m, xxx
% called by xxx
%

clc
clear
close all
format short, format compact

%-----------------------------------

% add path to additional matlab scripts
path(path,[pwd '/matlab_scripts']);

%-----------------------------------

% USER OUTPUT DIRECTORIES
dir0 = '/home/carltape/ADJOINT_TOMO/ADJOINT_TOMO_OUTPUT/';
pdir = [dir0 'model_plot_matlab_OUTPUT/'];
mdir = [dir0 'OUTPUT_subspace_save_m16/'];

% KEY: vertical or horizontal cross-sections
ivert = 0;
iplot_vert = 0;

ihorz = 0;
iplot_horz = 1;
iplot_horz_coverage = 0;
iplot_horz_ker = 0;

iwrite = 1;
ifigure = 0;
imask = 1;      % write out masked versions with NaN (=1) or not (=0)

% UTM zone
szone = '11S';

%-----------------------------------

% get a set of plotting gridpoints to extract a function value at a
% given x,y,z from a set of .bin files

% bounds for SPECFEM3D simulations
ax0 = [-121.6 -114.7 32.2 36.8];
ax0_utm = axes_utm2ll(ax0,szone,0);
%ax0_utm = [ 0.3271    0.4698    3.6792    3.7914] * 1e6;   % zoom
xmin = ax0_utm(1); xmax = ax0_utm(2);
ymin = ax0_utm(3); ymax = ax0_utm(4);

 % outline of SPECFEM region
[xbox_utm,ybox_utm] = boxpts(ax0_utm,500);
[xbox_ll,ybox_ll] = utm2ll(xbox_utm,ybox_utm,szone,1);

figure; plot(xbox_utm, ybox_utm, '.');

%-----------------------------------
% HORIZONTAL CROSS-SECTIONS

if ihorz == 1
    
    irun = 2;   % 1 for 1 km depth increment
                % 2 for 0.5 km depth increment
    dir1 = [pdir 'horz_' sprintf('%2.2i',irun) '/'];

    switch irun
        case 1
            depvec = 1e3 * [ 0:-1:-40];  % KEY
            numx = 300;         % key command for resolution
            xvec = linspace(xmin,xmax,numx);
            dx = xvec(2) - xvec(1);
            yvec = [ymin:dx:ymax];
            numy = length(yvec);
        case 2
            depvec = 1e3 * [ 0:-0.5:-40];  % KEY
            dx = 2000;         % key command for resolution
            xvec = [xmin:dx:xmax];
            yvec = [ymin:dx:ymax];
            numx = length(xvec);
            numy = length(yvec);
    end
    ndep = length(depvec);
    [X,Y] = meshgrid(xvec,yvec);
    Xvec = X(:);
    Yvec = Y(:);
    N = length(Xvec);
    
    % convert back to lat-lon
    disp('Entering utm2ll.m...');
    [Xvec_ll,Yvec_ll] = utm2ll(Xvec,Yvec,szone,1);
    figure; plot(Xvec_ll, Yvec_ll, '.','markersize',2);

    disp(sprintf('horizontal mesh is %i (x) by %i (y) for %i points total',numx,numy,N));
    disp(sprintf('horizontal spacing of plotting points is %.2f meters',dx));

    %----------------------------

    if iwrite==1
        % list of depths for each index
        fid = fopen([dir1 'horz_xc_cuts'],'w');
        for ii=1:ndep
            fprintf(fid,'%10i%16i\n',ii,depvec(ii));
        end
        fclose(fid);

        % reference horizontal surface
        fid = fopen([dir1 'horz_xc_lonlat'],'w');
        for kk=1:N
            fprintf(fid,'%16.6e%16.6e%16.6e%16.6e\n',Xvec_ll(kk),Yvec_ll(kk),Xvec(kk),Yvec(kk));
        end
        fclose(fid);

        % xyz files
        for ii=1:ndep
            zdep = depvec(ii)
            filename = sprintf('horz_xc_%3.3i',ii);
            disp(['writing file: ' dir1 filename]);
            fid = fopen([dir1 filename],'w');
            for kk=1:N
                fprintf(fid,'%16.6e%16.6e%16.6e\n',Xvec(kk),Yvec(kk),zdep);
            end
            fclose(fid);
        end
    end
end

%----------------------------
% VERTICAL CROSS-SECTIONS

if ivert ==1

    % read in sources and receivers
    recfile = '/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem';
    [rlon,rlat,relev,rburial,stnm,netwk] = read_station_SPECFEM(recfile);
    stanet = merge_arrays(stnm,netwk,'.');
    srcfile = '/home/carltape/results/SOURCES/socal_16/EIDs_lonlat_otime';
    [~,eid,~,~,sdep,~,~,slon,slat] = textread(srcfile,'%f%f%s%s%f%f%f%f%f');

    % irun
    % 1 src-rec
    % 2 SAF + 3
    % 3 Garlock
    % 4 LON
    % 5 LAT
    irun = 5;   % KEY: which set of cross-sections to make
    stirun = sprintf('%2.2i',irun);
    dir1 = [pdir 'vert_' stirun '/'];   % KEY COMMAND
    
    switch irun
        
        case 1      % 1: PICK SOURCE-RECEIVER PAIRS
            % specify the src-rec pairs and find their lon,lat,dep
            eids = [9983429 14096196 9968977 14179736 14383980 9818433 14418600 9967901 ...
                14383980 9703873 14236768 10992159 14138080 14186612 14186612 9688709 ...
                14077668 14095540 9627721 14236768 10215753 10215753 14096196 14096196 ...
                9983429 10006857 9173365 14095540 14077668 14095628 9674049 3321590 ...
                3320736 14179292 9674049 14418600 9173365 9155518 9718013 14179736 ...
                9154092 14169456 9045109 14418600 3320736 13935988 14179292 14236768 ...
                14383980 9695397 9173365 9882325 14077668 10006857 9967901 9695397 ...
                14155260 14077668 14096196 9818433 10006857 10215753 14151344 9173365 ...
                13935988 14155260 10097009 9882329 14169456 9674213 9753485 9096972 ...
                14000376 10370141 14383980 9093975 12659440 14077668 9753485 14155260 ...
                14155260 14155260 14155260 14155260 14155260 14155260];
            recs = {'LDF.CI','WGR.CI','HEC.CI','LAF.CI','STC.CI','CLC.CI','HEC.CI','OSI.CI',...
                'SMS.CI','RVR.CI','EDW2.CI','SMM.CI','BZN.AZ','PSR.CI','FMP.CI','BRE.CI',...
                'LDF.CI','LCG.CI','FMP.CI','WSS.CI','LCG.CI','HAST.TA','SMM.CI','SCI2.CI',...
                'DAN.CI','LGB.CI','FRD.AZ','SNCC.CI','SNCC.CI','SNCC.CI','SNCC.CI','SNCC.CI',...
                'SNCC.CI','SNCC.CI','LDF.CI','LDF.CI','LDF.CI','LDF.CI','LDF.CI','LDF.CI',...
                'LDF.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI',...
                'ISA.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI','ISA.CI','CWC.CI',...
                'SCZ2.CI','SCZ2.CI','SCZ2.CI','SCZ.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI',...
                'LGU.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI',...
                'LGU.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI','LGU.CI','DPP.CI','TEH.CI',...
                'LDR.CI','LRL.CI','CLC.CI','SLA.CI','LMR2.CI','DSC.CI'};
            
            slons = zeros(length(eids),1);
            slats = zeros(length(eids),1);
            sdeps = zeros(length(eids),1);
            nrays = length(eids);
            for ii=1:nrays
                imatch = find(eids(ii) == eid);
                slons(ii) = slon(imatch);
                slats(ii) = slat(imatch);
                sdeps(ii) = sdep(imatch);
            end
            rlons = zeros(nrays,1); rlats = zeros(nrays,1); rdeps = zeros(nrays,1);
            for ii=1:nrays
                imatch = find(strcmp(recs(ii),stanet)==1);
                rlons(ii) = rlon(imatch);
                rlats(ii) = rlat(imatch);
            end

            % specify parameters for sampling
            preseg_vec = 100*ones(nrays,1);
            postseg_vec = 100*ones(nrays,1);
            
        case 2            % 2: San Andreas
            
            % select profiles perpendicular to the SAF -- saf.m
            saflonlat_cuts = [
                -120.9420   36.3617 -118.4337   38.1771
                -120.7410   36.1754 -118.2194   37.9734
                -120.5640   36.0081 -117.9244   37.6869
                -120.3880   35.8542 -117.9540   37.7197
                -120.2380   35.6755 -117.8501   37.5748
                -119.9610   35.4046 -117.3811   37.1230
                -119.7050   35.1649 -117.3195   37.0532
                -119.5490   35.0317 -117.3772   37.0805
                -119.4150   34.9367 -117.4435   37.1133
                -119.2290   34.8725 -118.1995   37.4393
                -119.0590   34.8386 -117.5635   37.2505
                -118.9380   34.8158 -118.0720   37.4215
                -118.7480   34.7690 -117.7609   37.3463
                -118.5930   34.7257 -117.2598   37.1985
                -118.4080   34.6632 -117.0520   37.1275
                -118.1260   34.5614 -116.3847   36.8594
                -117.8940   34.4664 -116.2315   36.8018
                -117.6070   34.3428 -115.8980   36.6545
                -117.4670   34.2764 -115.7606   36.5887
                -117.2620   34.1669 -115.5181   36.4592
                -116.9810   34.0651 -115.5469   36.4959
                -116.7950   33.9496 -116.7731   36.6475
                -116.5660   33.9238 -115.3337   36.4268
                -116.3040   33.8099 -114.4192   36.0215
                -116.0630   33.6362 -113.9734   35.7163
                -115.8780   33.4662 -113.7552   35.5201
                -115.6280   33.2626 -113.3162   35.1667
                -115.5940   33.0397 -113.1488   34.8184
                -115.4120   32.7096 -112.9413   34.4540 ];
            nsafcuts = length(saflonlat_cuts);
            
            slons = [-118.1 -118.1 -118.1 saflonlat_cuts(:,1)'];
            slats = [33.9 33.9 33.9 saflonlat_cuts(:,2)'];
            rlons = [-117 -119.4 -120 saflonlat_cuts(:,3)'];
            rlats = [35 35.4 35 saflonlat_cuts(:,4)'];
            nrays = length(slons);
            preseg_vec = [100 150 150 200*ones(1,nsafcuts)];
            postseg_vec = [100 150 150 0*ones(1,nsafcuts)];
            sdeps = zeros(nrays,1);
            
        case 3            % 3: Garlock
            
            % select profiles perpendicular to the Garlock -- saf.m
            garlonlat_cuts = [
                 -118.7309   34.8803 -117.4704   33.4150
                 -118.5238   34.9737 -117.6324   33.3333
                 -118.1984   35.1372 -116.8821   33.7032
                 -118.0406   35.2404 -116.5687   33.9113
                 -117.8138   35.3766 -116.5167   33.9286
                 -117.4834   35.4934 -116.6783   33.8211
                 -117.1826   35.5674 -116.6713   33.8186
                 -116.9706   35.6063 -116.9243   33.8081
                 -116.6402   35.5966 -116.6391   33.7980
                 -116.4528   35.5966 -116.2956   33.8026  ];
            ngarcuts = length(garlonlat_cuts);
            
            slons = garlonlat_cuts(:,1)';
            slats = garlonlat_cuts(:,2)';
            rlons = garlonlat_cuts(:,3)';
            rlats = garlonlat_cuts(:,4)';
            nrays = length(slons);
            preseg_vec = 200*ones(1,ngarcuts);
            postseg_vec = 0*ones(1,ngarcuts);
            sdeps = zeros(nrays,1);
            
        case 4
            % cuts of constant longitude
            lon_cuts = [
             -115.0860   32.2367 -115.0860   37.1830
             -115.5080   32.9593 -115.5080   37.9056
             -115.9960   33.5937 -115.9960   38.5400
             -116.4850   33.8932 -116.4850   38.8395
             -116.9810   34.0651 -116.9810   39.0114
             -117.5030   34.2957 -117.5030   39.2420
             -118.0110   34.5118 -118.0110   39.4581
             -118.5090   34.6976 -118.5090   39.6439
             -118.9980   34.8209 -118.9980   39.7672
             -119.4850   34.9868 -119.4850   39.9331
             -119.9980   35.4483 -119.9980   40.3946
             -120.4920   35.9436 -120.4920   40.8899
             -120.9890   36.4062 -120.9890   41.3525   ];
            nloncuts = length(lon_cuts);

            slons = lon_cuts(:,1)';
            slats = lon_cuts(:,2)';
            rlons = lon_cuts(:,3)';
            rlats = lon_cuts(:,4)';
            nrays = length(slons);
            preseg_vec = 250*ones(1,nloncuts);
            postseg_vec = 0*ones(1,nloncuts);
            sdeps = zeros(nrays,1);
            
        case 5
            % cuts of constant latitude
            lat_cuts = [
                 -115.2100   32.5132 -114.7000   32.5132
                 -115.5940   33.0397 -114.7000   33.0397
                 -115.9090   33.4981 -114.7000   33.4981
                 -116.8660   34.0004 -114.7000   34.0004
                 -117.9800   34.5007 -114.7000   34.5007
                 -119.4850   34.9868 -114.7000   34.9868
                 -120.0410   35.4953 -114.7000   35.4953
                 -120.5640   36.0081 -114.7000   36.0081
                 -121.0880   36.5124 -114.7000   36.5124   ];
            nlatcuts = length(lat_cuts);

            slons = lat_cuts(:,1)';
            slats = lat_cuts(:,2)';
            rlons = lat_cuts(:,3)';
            rlats = lat_cuts(:,4)';
            nrays = length(slons);
            preseg_vec = 300*ones(1,nlatcuts);
            postseg_vec = 0*ones(1,nlatcuts);
            sdeps = zeros(nrays,1);
    end
    
    fid = fopen([dir1 'ALL_rays_' stirun],'w');
    aziall = zeros(nrays,2);
    for ii = 1:nrays
        azi1 = azimuth(rlats(ii),rlons(ii),slats(ii),slons(ii));
        azi2 = azimuth(slats(ii),slons(ii),rlats(ii),rlons(ii));
        aziall(ii,:) = [azi1 azi2];
        if irun==1
            fprintf(fid,'%3i%12i%10s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',...
                ii,eids(ii),recs{ii},slons(ii),slats(ii),rlons(ii),rlats(ii),aziall(ii,:));
        else
            fprintf(fid,'%3i%12i%10s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',...
                ii,NaN,'NAN',slons(ii),slats(ii),rlons(ii),rlats(ii),aziall(ii,:));            
        end
    end
    fclose(fid);
    
    dx = 1;     % horizontal sampling interval, km
    zmax = 5*1e3; zmin = -45*1e3; znum = 100;
   
%     % update the modified file if you add any columns (i.e., faults)
%     if 0==1
%         dir1 = '/home/carltape/results/PLOTTING/SAVE/vert_01/';
%         [i1,i2,i3,i4,i5,f1,f2,f3,f4,f5,f6] = ...
%             textread([dir1 'ALL_rays_fault_positions_mod_old'],'%f%f%s%f%f%f%f%f%f%f%f');
%         [i1b,i2b,i3b,i4b,i5b,f1b,f2b,f3b,f4b,f5b,f6b,f7b,f8b,f9b] = ...
%             textread([dir1 'ALL_rays_fault_positions'],'%f%f%s%f%f%f%f%f%f%f%f%f%f%f');
%         fid = fopen([dir1 'ALL_rays_fault_positions_new'],'w');
%         stfmt = repmat('%10.0f',1,9);
%         for ii = 1:length(i1b)
%             fprintf(fid,['%3i%12i%10s%10.1f%10.1f' stfmt '\n'],...
%                 i1b(ii),i2b(ii),i3b{ii},i4b(ii),i5b(ii),...
%                 f1(ii),f2(ii),f3(ii),f4(ii),f5(ii),f6(ii),...
%                 f7b(ii),f8b(ii),f9b(ii));
%         end
%         fclose(fid);
%     end
    
    % load six faults
    ffile = '/home/carltape/gmt/faults/socal_faults_xc.lonlat';
    [data,ilabs,labels,finds] = read_delimited_file(ffile,2);
    nfault = length(finds);
    DMIN_KM = 10;    % min distance must be DMIN_KM from the fault to mark the label
    
    % compute distances to six selected faults
    fid0 = fopen([dir1 'ALL_rays_' stirun '_fault_positions'],'w');
    
    imin = 1; imax = nrays;
    %imin = 80; imax = nrays;
    %imin = nrays; imax = imin;
    
    for ii = imin:imax
        
        % source and receiver points
        lon1 = slons(ii); lat1 = slats(ii);
        lon2 = rlons(ii); lat2 = rlats(ii);
        
        if irun==1
            stit = sprintf('event %i to %s',eids(ii),recs{ii});
        else
            stit = 'model cut';
        end
        disp(sprintf('%s -- %i out of %i to %i',stit,ii,imin,imax));
        
%         switch ii
%             case 1
%                 % this example has input points going outside the mesh:
%                 % the bottom, top, and south side
%                 znum = 100;
%                 hnum = 200;
%                 zmax = 5*1e3; zmin = -70*1e3;
%                 lon1 = -119.6; lat1 = 35.6; lon2 = -115.9; lat2 = 31.2;
%                 %rng1 = 600; az1 = 144;
%             case 2
%                 znum = 100;
%                 hnum = 100;
%                 zmax = 5*1e3; zmin = -45*1e3;
%                 lon1 = -119.6; lat1 = 34.6; lon2 = -117.4; lat2 = 33.4;
%                 %rng1 = 220; az1 = 110;
%         end

        % depth sampling
        depvec = linspace(zmin,zmax,znum)';

        % get points on the ray path
        preseg = preseg_vec(ii);
        postseg = postseg_vec(ii);
        [lon_gc1,lat_gc1,rdist,azis] = raypath(lon1,lat1,lon2,lat2,preseg,postseg,dx);  % KEY COMMAND
        [lon_gc1_utm,lat_gc1_utm] = utm2ll(lon_gc1,lat_gc1,szone,0);
        
        % exclude the points on the ray that are outside the mesh,
        % then convert back to lat-lon
        inds = getsubset(lon_gc1_utm,lat_gc1_utm,ax0_utm);
        rdist = rdist(inds);
        azis = azis(inds);
        lon_gc1_utm = lon_gc1_utm(inds);
        lat_gc1_utm = lat_gc1_utm(inds);
        [lon_gc1,lat_gc1] = utm2ll(lon_gc1_utm,lat_gc1_utm,szone,1);
        hnum = length(lon_gc1);
        
        % compute distance to each of the six faults
        dmins = 1e9*ones(nfault,1);
        rinds = zeros(nfault,1);
        for kk = 1:nfault
            indf = [finds(kk,1):finds(kk,2)];
            flon = data(indf,1); flat = data(indf,2); nf = length(flon);
        
            dminall = 1e6*ones(nf,1);
            iminall = zeros(nf,1);
            for ij = 1:nf
                [dtemp,itemp] = min(deg2km(distance(...
                    flon(ij)*ones(hnum,1),flat(ij)*ones(hnum,1),lon_gc1,lat_gc1)));
                dminall(ij) = dtemp;
                iminall(ij) = itemp;
            end
            [gmin,gind] = min(dminall);
            rinds(kk) = iminall(gind);
            dmins(kk) = dminall(gind);
        end
        
        rcuts = 9999*ones(1,nfault);
        igood = find(dmins < DMIN_KM);
        rcuts(igood) = rdist(rinds(igood));
        
        % write the fault cuts to file
        stfmt = repmat('%10.0f',1,nfault);
        if irun==1
            fprintf(fid0,['%3i%12i%10s%10.1f%10.1f' stfmt '\n'],...
                ii,eids(ii),recs{ii},aziall(ii,:),rcuts);
        else
            fprintf(fid0,['%3i%12i%10s%10.1f%10.1f' stfmt '\n'],...
                ii,NaN,'NaN',aziall(ii,:),rcuts);

        end
        
%         % get points on the ray path
%         %[lat_gc1,lon_gc1] = track1(lat1,lon1,az1,km2deg(rng1),[],'degrees',hnum);
%         [lat_gc1,lon_gc1] = gcwaypts(lat1,lon1,lat2,lon2,hnum-1);
%         [lon_gc1_utm,lat_gc1_utm] = utm2ll(lon_gc1,lat_gc1,szone,0);
% 
%         % compute the distance in km from one point to the end --
%         % this is used in the 2D plots
%         [rng,az] = distance(lat1*ones(hnum,1),lon1*ones(hnum,1),lat_gc1,lon_gc1);
%         rdist = deg2km(rng);
        
        figure; hold on;
        plot(xbox_ll, ybox_ll, '.');
        plot(lon_gc1,lat_gc1,'b.');
        plot(lon_gc1(rinds),lat_gc1(rinds),'ko','markersize',10);
        for kk=1:nfault, indf = [finds(kk,1):finds(kk,2)]; plot(data(indf,1),data(indf,2)); end
        text(lon_gc1(1),lat_gc1(1),num2str(1),'fontsize',16);
        text(lon_gc1(hnum),lat_gc1(hnum),num2str(hnum),'fontsize',16);
        plot(lon1,lat1,'rp',lon2,lat2,'ko','markersize',20);
        title(stit);

%         % reference points at the surface
%         % COLUMNS: dist from start (km), utm-x, utm-y, lon, lat
%         filename = sprintf('socal_xc_%3.3i_ref',ii);
%         fid = fopen([pdir filename],'w');
%         for kk=1:hnum    % loop over horizontal points in the ray path
%         	fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e\n',rdist(kk),...
%                 lon_gc1_utm(kk),lat_gc1_utm(kk),lon_gc1(kk),lat_gc1(kk));
%         end
%         fclose(fid);
        
        if iwrite == 1
            % xyz files
            filename = sprintf('vert_%2.2i_xc_%3.3i_input',irun,ii);
            disp(['writing file: ' dir1 filename]);
            fid = fopen([dir1 filename],'w');
            for kk=1:hnum    % loop over horizontal points in the ray path
                for jj=1:znum
                    fprintf(fid,'%16.6e%16.6e%16.6e%16.6e\n',...
                        lon_gc1_utm(kk),lat_gc1_utm(kk),depvec(jj),1000*rdist(kk));
                end
            end
            fclose(fid);

            % ray path
            filename = sprintf('vert_%2.2i_xc_%3.3i_ray_path',irun,ii);
            fid = fopen([dir1 filename],'w');
            for kk=1:hnum    % loop over horizontal points in the ray path
                fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e\n',...
                    lon_gc1(kk),lat_gc1(kk),lon_gc1_utm(kk),lat_gc1_utm(kk),1000*rdist(kk));
            end
            fclose(fid);
            
            % lon-lat points of source and receiver
            filename = sprintf('vert_%2.2i_xc_%3.3i_A_B',irun,ii);
            fid = fopen([dir1 filename],'w');
            [dq,azq] = distance(slats(ii),slons(ii),rlats(ii),rlons(ii));   % distance to receiver
            iflip = 0;
            %if azq > 180, iflip = 1; end
            if rlons(ii) < slons(ii), iflip = 1; end
            fprintf(fid,[repmat('%16.6e',1,9) '%4i\n'],...
                slons(ii),slats(ii),0,sdeps(ii),rlons(ii),rlats(ii),deg2km(dq),0,azq,iflip);
            fclose(fid);
        end
        
    end   % for ii=1:nrays
    fclose(fid0);       % fault cut positions
    
end

%----------------------------
% PLOTTING CROSS-SECTIONS

if iplot_vert==1
    colors;
    
    %irun = 1; pmax = 86;
    %irun = 2; pmax = 32;
    %irun = 3; pmax = 10;
    %irun = 4; pmax = 13;
    irun = 5; pmax = 9;
    
    stirun = sprintf('%2.2i',irun);
    dir1 = [pdir 'vert_' stirun '/'];
    smodel1 = 'm00';
    smodel2 = 'm16';
    modlab1 = 'vs';
    modlab2 = 'vb';
    smodeltag = [smodel2 '_' smodel1];

    pmin = 1;   % default
    %pmin = 80; pmax = 60;
    %pmin = 5; pmax = pmin;
    
    for ip = pmin:pmax
        stip = sprintf('%3.3i',ip);
        disp(sprintf('model %i out of %i',ip,pmax));
    
        % four possible combinations
        tag3 = [modlab1 '_' smodel1];
        tag4 = [modlab1 '_' smodel2];
        tag5 = [modlab2 '_' smodel1];
        tag6 = [modlab2 '_' smodel2];
        
        % load files needed for depth cross-sections
        fz0 = [dir1 'vert_' stirun '_xc_'];
        file0 = [fz0 stip '_A_B'];
        file1 = [fz0 stip '_ray_path'];
        file2 = [fz0 stip '_input'];
        file3 = [fz0 tag3 '_' stip '.gmt'];
        file4 = [fz0 tag4 '_' stip '.gmt'];
        file5 = [fz0 tag5 '_' stip '.gmt'];
        file6 = [fz0 tag6 '_' stip '.gmt'];

        [lon1,lat1,d1,z1,lon2,lat2,d2,z2,az,iflip] = textread(file0,'%f%f%f%f%f%f%f%f%f%f');
        [lon_gc1,lat_gc1,lon_gc1_utm,lat_gc1_utm,rdist] = textread(file1,'%f%f%f%f%f');
        [utm_x,utm_y,z,d] = textread(file2,'%f%f%f%f');
        [~,~,~,vs1,cdist] = textread(file3,'%f%f%f%f%f');
        [~,~,~,vs2,cdist] = textread(file4,'%f%f%f%f%f');
        [~,~,~,vb1,cdist] = textread(file5,'%f%f%f%f%f');
        [~,~,~,vb2,cdist] = textread(file6,'%f%f%f%f%f');
        n = length(vs1);
        
        %lon1 = lon_gc1(1); lat1 = lat_gc1(1);
        %lon2 = lon_gc1(end); lat2 = lat_gc1(end);
        
        % turn negative points into NaN values
        vs1(vs1 < 0) = NaN;
        vs2(vs2 < 0) = NaN;
        vb1(vb1 < 0) = NaN;
        vb2(vb2 < 0) = NaN;
        
        % compute Vp
        vp1 = zeros(n,1); vp2 = zeros(n,1);
        vp1 = sqrt( 4/3 * vs1.^2 + vb1.^2 );
        vp2 = sqrt( 4/3 * vs2.^2 + vb2.^2 );
        
        % compute the difference -- GMT does not seem to like operating on
        % NaN vectors, so we have to do it here and then write it out
        vs3 = log(vs2 ./ vs1);
        vb3 = log(vb2 ./ vb1);
        vp3 = log(vp2 ./ vp1);
        
        if ifigure == 1
            for kk = 1:3
                switch kk
                    case 1, cval = vs1;
                    case 2, cval = vs2;
                    case 3, cval = vs3;
                end

                % construct mesh for Matlab plotting
                xvec = d*1e-3;
                yvec = z*1e-3;
                zvec = cval;
                xlin  = linspace(min(xvec), max(xvec), 200);
                ylin  = linspace(min(yvec), max(yvec), 100);
                [X,Y] = meshgrid(xlin,ylin);

                % determine interpolated function using xvec,yvec input
                %whos xvec yvec zvec
                Z = griddata(xvec, yvec, zvec, X, Y, 'nearest');

                figure; nr=3; nc=1;
                %clims = [2000 8500];

                subplot(nr,nc,1); hold on;
                plot(xbox_ll, ybox_ll, '.');
                plot(lon_gc1,lat_gc1,'b.');
                plot(lon1,lat1,'rp','markersize',20);
                axis equal, axis tight

                subplot(nr,nc,2);
                pcolor(X,Y,Z); shading flat;
                xlabel('Distance along cross-section, km');
                ylabel('Depth, km');
                %caxis(clims);
                colorbar;

                subplot(nr,nc,3);
                pcolor(X,Y,Z); shading flat;
                xlabel('Distance along cross-section, km');
                ylabel('Depth, km');
                %caxis(clims);
                colorbar; axis equal, axis tight

                colormap(seis);
            end
        end
        
        if iwrite == 1
%             % write out a file with five columns: xdist, zdep, m00 (m/s), m16 (m/s), ln(m16/m00)
%             fid = fopen([fz0 modlab '_' smodeltag '_' stip '.dat'],'w');
%             for kk=1:length(d)
%                 fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e\n',...
%                     d(kk)*1e-3,z(kk)*1e-3,cval1(kk)*1e-3,cval2(kk)*1e-3,cval3(kk));
%             end
%             fclose(fid);
            
            % write out a file with 11 columns: xdist, zdep, m00 (m/s), m16 (m/s), ln(m16/m00)
            fid = fopen([fz0 'vsvbvp_' smodeltag '_' stip '.dat'],'w');
            for kk=1:length(d)
                fprintf(fid,[repmat('%16.6e',1,11) '\n'],...
                    d(kk)*1e-3,z(kk)*1e-3,vs1(kk)*1e-3,vs2(kk)*1e-3,vs3(kk),...
                                          vb1(kk)*1e-3,vb2(kk)*1e-3,vb3(kk),...
                                          vp1(kk)*1e-3,vp2(kk)*1e-3,vp3(kk) );
            end
            fclose(fid);
            %whos vs1 vs2 vs3 vb1 vb2 vb3 vp1 vp2 vp3, length(d)
        end
        
    end
    
end

%----------------------------
% PLOTTING MODEL CROSS-SECTIONS

if iplot_horz==1
    colors;
    
    irun = 2;       % index into set of horizontal cross sections
    stirun = sprintf('%2.2i',irun);
    dir1 = [pdir 'horz_' stirun '/'];
    smodel1 = 'm00';
    smodel2 = 'm16';    % m16, m01, etc
    modlab = 'vs';
    smodeltag = [smodel2 '_' smodel1];
    
    % KEY: TWO DIFFERENT KINDS OF SUMMED KERNELS
    % load number of window picks per event
    nwins_all = load([mdir 'm16/window_tomo/m16_window_picks_nwin']);
    nevent = length(nwins_all);
    nwin = sum(nwins_all);
    cfac = log10(1e-16);
    %covertag = 'coverage_sum_abs'; Nfac = nwin;
    covertag = 'coverage_sum_abs_nwin'; Nfac = nevent;
    cnorm = 10^cfac;
    
    % load files needed for horizontal cross-sections
    fname = [dir1 'horz_xc_cuts']; disp(['file : ' fname]);
    [inds,depvec] = textread(fname,'%f%f');
    fname = [dir1 'horz_xc_lonlat']; disp(['file : ' fname]);
    [xlon,ylat,xutm,yutm] = textread(fname,'%f%f%f%f');
    
    % for some reason matlab interpolation needs to work with smaller numbers,
    % so we subtract a reference value
    xutm0 = min(xutm);
    yutm0 = min(yutm);
    xutm = xutm - xutm0;
    yutm = yutm - yutm0;
    
    % replace values far west in regions of no coverage with NaN
    %xvpts = [-119 -122 -122 -119];
    %yvpts = [32 36 32 32];
    %imask = inpolygon(xlon,ylat,xvpts,yvpts);
          
    %figure; hold on;
    %plot(xlon,ylat,'b.');
    %plot(xlon(indmask),ylat(indmask),'ro');
    
    %pmin = 1; pmax = length(depvec);
    pmin = 21; pmax = pmin;
    %pmin = 1; pmax = 12;
    
    for ip = pmin:pmax
        stip = sprintf('%3.3i',ip);
        disp(sprintf('%i out of the set %i to %i',ip,pmin,pmax));

        % load output from cluster code
        tag1 = [modlab '_' smodel1];
        tag2 = [modlab '_' smodel2];
        file1 = [dir1 tag1 '/horz_xc_' tag1 '_' stip '.gmt'];
        disp(['file :' file1]);
        [~,~,~,cval1,cdist] = textread(file1,'%f%f%f%f%f');
        file2 = [dir1 tag2 '/horz_xc_' tag2 '_' stip '.gmt'];
        disp(['file :' file2]);
        [~,~,~,cval2,cdist] = textread(file2,'%f%f%f%f%f');
        
        if imask==1
            % load coverage mask
            tagm = [modlab '_m16_' covertag];
            %stip0=stip;
            %if ip==3, stip0 = '002'; end
            filem = [dir1 tagm '/horz_xc_' tagm '_' stip '.gmt'];
            [~,~,~,cvalm,cdist] = textread(filem,'%f%f%f%f%f');
        end

        % remove points above the topography -- these are set as negative
        % numbers in the interpolation algorithm
        cval1(cval1 < 0) = NaN;
        cval2(cval2 < 0) = NaN;
        
        % compute the difference -- GMT does not seem to like operating on
        % NaN vectors, so we have to do it here and then write it out
        cval3 = log(cval2 ./ cval1);
        
        %----------------------------------------------------
        % COVERAGE MASK
        
        if imask==1
            % take log10 of the field
            cvalmlog = log10(cvalm / Nfac);

            % replace points outside coverage with a mask
            indmask = find(cvalmlog <= cfac);
            cval1(indmask) = NaN;
            cval2(indmask) = NaN;
            cval3(indmask) = NaN;
        end
        
        %----------------------------------------------------     
        
        clabs = {smodel1,smodel2,['ln(' smodel2 '/' smodel1 ')']};
        for kk = 1:3
        %for kk = 2:2
            switch kk
                case 1, cval = cval1;
                case 2, cval = cval2;
                case 3, cval = cval3;
            end
            
            if ifigure ==1
                % construct mesh for plotting
                xlin  = linspace(min(xutm), max(xutm), 200);
                ylin  = linspace(min(yutm), max(yutm), 200);
                [X,Y] = meshgrid(xlin,ylin);

                % determine interpolated function using xvec,yvec input
                Z = griddata(xutm, yutm, cval, X, Y, 'nearest');

                figure;
                clims = [-1 1]*0.2;
                pcolor((X+xutm0)*1e-3,(Y+yutm0)*1e-3,Z);
                shading flat;
                %caxis(clims);
                xlabel('X distance, km');
                ylabel('Y distance, km');
                colorbar; axis equal, axis tight
                title(sprintf('Depth = %.1f km -- %s',depvec(ip)*1e-3,clabs{kk}));
                %title(sprintf('Depth = %.1f km (Vref = %.1f km/s)',depvec(ip)*1e-3,cnorm*1e-3));
                colormap(seis);
            end
        end

        if iwrite == 1
            % write out a file with five columns: lon, lat, m00 (m/s), m16 (m/s), ln(m16/m00)
            fname = [dir1 'horz_xc_' modlab '_' smodeltag '_' stip '_mask' num2str(imask) '.dat'];
            disp(['writing file : ' fname]);
            fid = fopen(fname,'w');
            for kk=1:length(xlon)
                fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e\n',...
                    xlon(kk),ylat(kk),cval1(kk),cval2(kk),cval3(kk));
            end
            fclose(fid);
        end
    end
    
end

%----------------------------
% PLOTTING COVERAGE CROSS-SECTIONS

if iplot_horz_coverage==1
    colors;
    
    irun = 2;
    stirun = sprintf('%2.2i',irun);
    dir1 = [pdir 'horz_' stirun '/'];
    modlab = 'vb_m16';
    %smodeltag = [smodel2 '_' smodel1];
    
    % load number of window picks per event
    nwins_all = load([mdir 'm16/window_tomo/m16_window_picks_nwin']);
    nevent = length(nwins_all);
    nwin = sum(nwins_all);
    
    % KEY: TWO DIFFERENT KINDS OF SUMMED KERNELS
    cfac = log10(1e-16);
    %smodeltag = 'coverage_sum_abs'; Nfac = nwin;
    smodeltag = 'coverage_sum_abs_nwin'; Nfac = nevent;
    cnorm = 10^cfac;
    
    % load files needed for depth cross-sections
    [inds,depvec] = textread([dir1 'horz_' stirun '_xc_cuts'],'%f%f');
    [xlon,ylat,xutm,yutm] = textread([dir1 'horz_' stirun '_xc_lonlat'],'%f%f%f%f');
    
    % replace values far west in regions of no coverage with NaN
    %xvpts = [-119 -122 -122 -119];
    %yvpts = [32 36 32 32];
    %indmask = inpolygon(xlon,ylat,xvpts,yvpts);
          
    %figure; hold on;
    %plot(xlon,ylat,'b.');
    %plot(xlon(indmask),ylat(indmask),'ro');
    
    pmin = 1; pmax = length(depvec);
    %pmin = 1; pmax = pmin;
    pmin = 13; pmax = length(depvec);
    
    for ip = pmin:pmax
        stip = sprintf('%3.3i',ip);
        disp(sprintf('%i out of %i - %i',ip,pmin,pmax));

        % load output from cluster code
        tag1 = [modlab '_' smodeltag];
        file1 = [dir1 tag1 '/horz_xc_' tag1 '_' stip '.gmt'];
        [~,~,~,cval1,cdist] = textread(file1,'%f%f%f%f%f');
        
        % for some reason matlab interpolation needs to work with smaller numbers,
        % so we subtract a reference value
        xutm0 = min(xutm);
        yutm0 = min(yutm);
        xutm = xutm - xutm0;
        yutm = yutm - yutm0;

        % remove points above the topography
        cval1(cval1 < 0) = NaN;
        
        % normalize by number of windows or number of events
        cval1 = cval1 / Nfac;
        
        % take -log10 of the field
        cval2 = log10(cval1);
        
        % replace points outside coverage with a mask
        indmask = find(cval2 <= cfac);
        cval3 = cval2; cval3(indmask) = NaN;
        cval1mask = cval1; cval1mask(indmask) = NaN;
        
        % convex hull -- convhull
        
        clabs = {tag1,['log10(' tag1 ')'],['log10(' tag1 ')']};
        for kk = 1:3
        %for kk = 2:2
            switch kk
                case 1, cval = cval1;
                case 2, cval = cval2;
                case 3, cval = cval3;
            end
            
            if ifigure ==1
                % replace negative values with NaN
                %ineg = find(cval < 0);
                %cval(ineg) = NaN;

                % construct mesh for plotting
                xlin  = linspace(min(xutm), max(xutm), 200);
                ylin  = linspace(min(yutm), max(yutm), 200);
                [X,Y] = meshgrid(xlin,ylin);

                % determine interpolated function using xvec,yvec input
                Z = griddata(xutm, yutm, cval, X, Y, 'nearest');

                figure;
                clims = [-1 1]*0.2;
                pcolor((X+xutm0)*1e-3,(Y+yutm0)*1e-3,Z);
                shading flat;
                %caxis(clims);
                xlabel('X distance, km');
                ylabel('Y distance, km');
                colorbar; axis equal, axis tight
                title(sprintf('Depth = %.1f km -- %s',depvec(ip)*1e-3,clabs{kk}));
                %title(sprintf('Depth = %.1f km (Vref = %.1f km/s)',depvec(ip)*1e-3,cnorm*1e-3));
                colormap(seis);
            end
        end

        if iwrite == 1
            % write out a file with five columns: lon, lat, log10(Ksum), log10(Ksum)+mask
            filetag = ['horz_xc_' modlab '_' smodeltag '_' stip];
            fid = fopen([dir1 filetag '.dat'],'w');
            for kk=1:length(xlon)
                fprintf(fid,'%16.6e%16.6e%16.6e%16.6e\n',...
                    xlon(kk),ylat(kk),cval1(kk),cval1mask(kk));
            end
            fclose(fid);
            
            % threshold value for mask
            fid = fopen([dir1 filetag '_mask_value.dat'],'w');
            fprintf(fid,'%.6e\n',cnorm);
            fclose(fid);
        end
    end
    
end

%----------------------------
% PLOTTING KERNEL CROSS-SECTIONS

if iplot_horz_ker == 1
    colors;
    
    irun = 2;
    stirun = sprintf('%2.2i',irun) ;
    dir1 = [pdir 'horz_' stirun '/'];
    
    %fsmooth = 'h003km_v001km';
    fsmooth = 'h002km_v001km';
    %fsmooth = 'h000km_v000km';
    ftags = {'Kmu_rayleigh','Kmu_reflection_one','Kmu_reflection_two'};
    
    % load files needed for horizontal cross-sections
    fname = [dir1 'horz_' stirun '_xc_cuts'];
    disp(['file : ' fname]);
    [inds,depvec] = textread([dir1 'horz_' stirun '_xc_cuts'],'%f%f');
    [xlon,ylat,xutm,yutm] = textread([dir1 'horz_' stirun '_xc_lonlat'],'%f%f%f%f');
    
    kvec = [1 5 9 13 17 21 25 29];
    
    %pmin = 1; pmax = length(depvec);
    %pmin = 13; pmax = pmin;
    %pmin = 13; pmax = length(depvec);
    
    for ik = 1:length(kvec)
        ip = kvec(ik)
        stip = sprintf('%3.3i',ip);
        %disp(sprintf('%i out of the set %i to %i',ip,pmin,pmax));
        
        % load output from cluster code
        file1 = [dir1 'kernels/horz_xc_' ftags{1} '_' fsmooth '_' stip '.gmt'];
        [~,~,~,K1,~] = textread(file1,'%f%f%f%f%f');
        file2 = [dir1 'kernels/horz_xc_' ftags{2} '_' fsmooth '_' stip '.gmt'];
        [~,~,~,K2,~] = textread(file2,'%f%f%f%f%f');
        file3 = [dir1 'kernels/horz_xc_' ftags{3} '_' fsmooth '_' stip '.gmt'];
        [~,~,~,K3,~] = textread(file3,'%f%f%f%f%f');

        % for some reason matlab interpolation needs to work with smaller numbers,
        % so we subtract a reference value
        xutm0 = min(xutm);
        yutm0 = min(yutm);
        xutm = xutm - xutm0;
        yutm = yutm - yutm0;
        
        %----------------------------------------------------     
        
        clabs = {'rayleigh-1','rayleigh-2','rayleigh-3'};
        for kk = 1:3
        %for kk = 2:2
            switch kk
                case 1, cval = K1;
                case 2, cval = K2;
                case 3, cval = K3;
            end
            
            if ifigure ==1
                % replace negative values with NaN
                %ineg = find(cval < 0);
                %cval(ineg) = NaN;

                % construct mesh for plotting
                xlin  = linspace(min(xutm), max(xutm), 200);
                ylin  = linspace(min(yutm), max(yutm), 200);
                [X,Y] = meshgrid(xlin,ylin);

                % determine interpolated function using xvec,yvec input
                Z = griddata(xutm, yutm, cval, X, Y, 'nearest');

                figure;
                clims = [-1 1]*0.2;
                pcolor((X+xutm0)*1e-3,(Y+yutm0)*1e-3,Z);
                shading flat;
                %caxis(clims);
                xlabel('X distance, km');
                ylabel('Y distance, km');
                colorbar; axis equal, axis tight
                title(sprintf('Depth = %.1f km -- %s',depvec(ip)*1e-3,clabs{kk}));
                %title(sprintf('Depth = %.1f km (Vref = %.1f km/s)',depvec(ip)*1e-3,cnorm*1e-3));
                colormap(seis);
            end
        end

        if iwrite == 1
            % write out a file with five columns: lon, lat, K1, K2, K3
            fid = fopen([dir1 'horz_xc_ker_' fsmooth '_' stip '.dat'],'w');
            for kk=1:length(xlon)
                fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e\n',...
                    xlon(kk),ylat(kk),K1(kk),K2(kk),K3(kk));
            end
            fclose(fid);
        end
    end
    
end

%====================================================================
