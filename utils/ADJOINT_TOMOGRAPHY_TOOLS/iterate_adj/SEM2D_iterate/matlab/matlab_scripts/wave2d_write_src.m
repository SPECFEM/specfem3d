%
% function 
% Carl Tape, 25-Jan-2010
%
% 
%
% calls xxx
% called by xxx
%

function wave2d_write_src(filename,xs0,ys0,ts,xs,ys,ts_res,xs_res,ys_res)

disp(['writing to ' filename]);

% ! sources for synthetics
% write(20,'(8e20.10)') xtemp, ztemp, &
% m_src_syn_vec(itemp1), m_src_syn_vec(itemp2), m_src_syn_vec(itemp3), &
% (m_src_syn_vec(itemp1) - m_src_dat_vec(itemp1)), &
% (m_src_syn_vec(itemp2) - m_src_dat_vec(itemp2)), &
% m_src_syn_vec(itemp3) - m_src_dat_vec(itemp3)

fid = fopen(filename,'w');
for ii = 1:length(xs0)
    fprintf(fid,'%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e%20.10e\n',...
        xs0(ii), ys0(ii),...
        ts(ii), xs(ii), ys(ii),...
        ts_res(ii), xs_res(ii), ys_res(ii) );
end
fclose(fid);

%=========================================================