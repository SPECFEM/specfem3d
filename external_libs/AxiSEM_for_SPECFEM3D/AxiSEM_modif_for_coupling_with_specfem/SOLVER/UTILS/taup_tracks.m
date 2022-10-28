function run_taup(distance, phase, depth, timestep, src_shift)

cmd_pierce = sprintf('%d,', 10:1:2891);
cmd = sprintf('~/sources/TauP-2.1.2/bin/taup_pierce -nodiscon -h %d -ph %s -deg %d -pierce %s2891 >Pierce.txt', ...
              depth, phase, distance, cmd_pierce);
system(cmd)

% xc = 6371*sin(linspace(0,2*pi,1000));
% yc = 6371*cos(linspace(0,2*pi,1000));

% figure; hold on
% plot(xc, yc)

src_shift = fix(src_shift/timestep) * timestep;


%%
nheader = 1;
nphase  = 0;
% ival = 1;
fid_pierce = fopen('Pierce.txt');
while (true) 
    fgetl(fid_pierce)
    if (ans==-1)
        break
    end
    clear r theta timevals
    try 
%             [r(ival), theta(ival), timevals(ival)] = fscanf(fid_pierce, '%f ', 3);
            res = textscan(fid_pierce, '%f %f %f ');
    catch 
        break
    end
    r = res{2};
    theta = res{1};
    timevals = res{3} + src_shift;
    nlayer = length(r);
    
    theta(1:nlayer-1) = (theta(1:nlayer-1) + theta(2:nlayer)) / 2;

    
    nphase = nphase + 1;
    sprintf('Found phase number %d', nphase);
%     nheader = nheader + size(pierce.data,1) + 1;

    y = 1000*(6371-r) .* cosd(theta);
    x = 1000*(6371-r) .* sind(theta);
    
    traveltime = max(timevals);


    % plot(x,y)
    %%
    xdmf_start_file = ['<?xml version="1.0" ?>\n', ...
                      '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n',...
                      '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n',...
                      '<Domain>\n',...
                      '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">\n'];

    xdmf_start_time = ['   <Grid Name="%04d" GridType="Uniform">     \n', ...
                       '     <Time Value="%8.2f" />                \n', ...
                       '         <Topology TopologyType="Polyline"  \n',...
                       '          NodesPerElement="2" NumberOfElements="%d"> \n', ...
                       '             <DataItem Format="XML"         \n', ...
                       '                 Dimensions="%d 2"          \n', ...
                       '                 DataType="Int">            \n'];
    xdmf_data_topo  =  '                  %d %d \n';          
    xdmf_interm     = ['             </DataItem>                   \n', ...
                       '         </Topology>                         \n', ... 
                       '         <Geometry Type="XYZ">               \n', ...
                       '           <DataItem Format="XML" Dimensions="%d 3">\n'];
    xdmf_data       =  '             %d %d 0.0\n';
    xdmf_end_time   = ['	       </DataItem>                       \n', ...
                       '         </Geometry>                         \n', ...
                       '    </Grid>\n\n'];
    xdmf_end_file   =  '</Grid>\n</Domain>\n</Xdmf>\n';

    %%

    out = sprintf(xdmf_start_file);

    for itime = timestep:timestep:traveltime
        timeval = itime;
        np_at_time = find(timevals<timeval, 1, 'last');

        if (np_at_time>0)     
            out = [out, sprintf(xdmf_start_time, itime, timeval, np_at_time-2, np_at_time-1)];
            out = [out, sprintf(xdmf_data_topo,  [(1:np_at_time-1); (2:np_at_time)]) ];
            out = [out, sprintf(xdmf_interm,     np_at_time)];
            out = [out, sprintf(xdmf_data, [x(1:np_at_time)'; y(1:np_at_time)'])];
            out = [out, sprintf(xdmf_end_time)];
        end
    end
    out = [out, sprintf(xdmf_end_file)];
    fnam = sprintf('line_%s_%03d_%03d_%1d.xdmf', phase, distance, depth, nphase);
    fid = fopen(fnam, 'w');
    fprintf(fid, out);
    fclose(fid);


end
    
fclose(fid_pierce);
end
    
               
               