%
% misfit_gmt.m
% CARL TAPE, 03-Sept-2008
% printed xxx
%
%
% calls xxx
% called by xxx
%

%-------------------------

function norms = misfit_gmt(odir,Tmin,Tmax,smodel1,smodel2,eid,stnm,netwk)

disp('running misfit_gmt.m...');


%-------------------------

ifigure = 0;

Ttag = sprintf('T%3.3i_T%3.3i',Tmin,Tmax);
ssuffix1 = ['semd.sac.d.' Ttag '.' smodel1];
ssuffix2 = ['semd.sac.d.' Ttag '.' smodel2];
dsuffix = ['sac.d.' Ttag];

comps = {'Z','R','T'};
%chans = {'HL','BL','HN','BN','HH','BH'};
chans = {'HH','BH'};
ncomp = length(comps);
nchan = length(chans);

%dir0 = pwd;
%odir = [dir0 '/OUTPUT_MISFIT/EVENTS/' eid '/' stnm '.' netwk '/'];

%-------------------------

fhits = zeros(ncomp,3);
ffiles = repmat(cellstr(''),ncomp,3);
for icomp = 1:ncomp
    for ichan = 1:nchan
        
        comp = [chans{ichan} comps{icomp}];
        dfile = [odir eid '.' netwk '.' stnm '.' comp '.' dsuffix];
        sfile1 = [odir stnm '.' netwk '.' comp '.' ssuffix1];
        sfile2 = [odir stnm '.' netwk '.' comp '.' ssuffix2];
    
        if exist(dfile,'file')
            fhits(icomp,1) = fhits(icomp,1)+1;
            ffiles{icomp,1} = dfile;
        end
        if exist(sfile1,'file')
            fhits(icomp,2) = fhits(icomp,2)+1;
            ffiles{icomp,2} = sfile1;
        end
        if exist(sfile2,'file')
            fhits(icomp,3) = fhits(icomp,3)+1;
            ffiles{icomp,3} = sfile2;
        end
    
    end
end
fhits_sum = sum(fhits,2);

if ifigure==1
    figure; nr=5; nc=3;
    st1 = ['syn(' smodel1 ')'];
    st2 = ['syn(' smodel2 ')'];
end

norms = zeros(9,ncomp);

for icomp = 1:ncomp
%for icomp = 3   
    
    if fhits_sum(icomp) == 3
        dfile  = ffiles{icomp,1};
        sfile1 = ffiles{icomp,2};
        sfile2 = ffiles{icomp,3};

        [seisD,HdrData,tnu,pobj,timeD] = readsac(dfile,0,'l');
        [seisS1,HdrData,tnu,pobj,timeS1] = readsac(sfile1,0,'l');
        [seisS2,HdrData,tnu,pobj,timeS2] = readsac(sfile2,0,'l');

        ti = timeD(:);
        n = length(ti);

        % interpolate synthetics exactly onto the data time vector
        % (In theory, they should already be within two time-steps at this point.)
        s1 = interp1(timeS1,seisS1,ti,'linear','extrap');
        s2 = interp1(timeS2,seisS2,ti,'linear','extrap');

        % residuals for the two models; also perturbation to seismogram
        res1 = s1 - seisD;
        res2 = s2 - seisD;
        spert = s2 - s1;
        
        ymax = 1.05*max(abs([seisD' s1' s2' res1' res2']));
        tmin = min(ti); tmax = max(ti);
        ax0 = [tmin tmax -ymax ymax];

        % compute nine different norms -- IS THIS EXACTLY WHAT WE WANT?
        norms(1,icomp) = norm(seisD);
        norms(2,icomp) = norm(s1);
        norms(3,icomp) = norm(s2);
        norms(4,icomp) = norm(res1);
        norms(5,icomp) = norm(res2);
        norms(6,icomp) = norm(s2-s1);
        norms(7,icomp) = norm(res1)^2 / norm(seisD)^2;             % chi^2
        norms(8,icomp) = norm(res2)^2 / norm(seisD)^2;             % chi^2
        norms(9,icomp) = 100*(1 - (norm(res2)^2 / norm(res1)^2));  % VR
        
        if ifigure==1
            subplot(nr,nc,icomp);
            plot(ti,s1,'r',ti,seisD,'b'); legend(st1,'data'); axis(ax0);
            subplot(nr,nc,nc*1+icomp);
            plot(ti,res1,'k'); legend([st1 ' - data']); axis(ax0);
            subplot(nr,nc,nc*2+icomp);
            plot(ti,s2,'r',ti,seisD,'b'); legend(st2,'data'); axis(ax0);
            subplot(nr,nc,nc*3+icomp);
            plot(ti,res2,'k'); legend([st2 ' - data']); axis(ax0);
            subplot(nr,nc,nc*4+icomp);
            plot(ti,spert,'k'); legend([st2 ' - ' st1]); axis(ax0);
        end

        % WRITE FILES FOR GMT PLOTTING
        filename = [odir 'GMT_time_series_' Ttag '_' comps{icomp} '.dat'];
        fid = fopen(filename,'w');
        for ii = 1:n
            fprintf(fid,[repmat('%12.4e',1,7) '\n'],...
                ti(ii),seisD(ii),s1(ii),s2(ii),res1(ii),res2(ii),spert(ii) );   
        end
        fclose(fid);
        
        filename = [odir 'GMT_axes_' Ttag '_' comps{icomp} '.dat'];
        fid = fopen(filename,'w');
            fprintf(fid,[repmat('%12.4e',1,4) '\n'],...
                ax0(1),ax0(2),ax0(3),ax0(4) );   
        fclose(fid);
        
        filename = [odir 'GMT_norms_' Ttag '_' comps{icomp} '.dat'];
        fid = fopen(filename,'w');
        for kk = 1:9, fprintf(fid,'%12.4e\n',norms(kk,icomp) ); end
        fclose(fid);
        
    end
    
end

if ifigure == 1 
    orient tall, wysiwyg, fontsize(6)
end
    
return

%-----------------------------------
% if all three components exist, then compute additional time series for plots

if sum(fhits_sum) == 9
    
    % compute the time vector to use for all components
    mint = zeros(ncomp,3);
    maxt = zeros(ncomp,3);
    for icomp = 1:ncomp
        dfile  = ffiles{icomp,1};
        sfile1 = ffiles{icomp,2};
        sfile2 = ffiles{icomp,3};
        [seisD,HdrData,tnu,pobj,timeD] = readsac(dfile,0,'l');
        [seisS1,HdrData,tnu,pobj,timeS1] = readsac(sfile1,0,'l');
        [seisS2,HdrData,tnu,pobj,timeS2] = readsac(sfile2,0,'l');
        
        mint(icomp,1) = min(timeD);
        mint(icomp,2) = min(timeS1);
        mint(icomp,3) = min(timeS2);
        maxt(icomp,1) = max(timeD);
        maxt(icomp,2) = max(timeS1);
        maxt(icomp,3) = max(timeS2);
    end
    
    tmin = max(mint(:));
    tmax = min(maxt(:));
    ti = [tmin:0.05:tmax]';

end

%=======================================================================
