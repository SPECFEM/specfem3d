%
% function wave2d_plot_hessian(Hstr,Hsrc,H)
% Carl Tape, 16-Feb-2010
%
% Plot the Hessian for the source subspace projection method.
%
% calls xxx
% called by wave2d_subspace.m, wave2d_cg_run.m
%

function wave2d_plot_hessian(Hstr,Hsrc,H)

disp(' properties of Hessian (min, median, mean(abs), max, std):');
stH = sprintf('min %.2e, median %.2e, mean(abs) %.2e, max %.2e, std %.2e',...
            min(H(:)), median(H(:)), mean(abs(H(:))), max(H(:)), std(H(:)));
disp(stH);

%if ~and(INV_STRUCT==1, INV_SOURCE==1)
if ~and(norm(Hstr)==0, norm(Hsrc)==0)
    
    figure;
    %pcolor(H); shading flat;
    imagesc(H);
    xlabel('Source index'); ylabel('Source index');
    title({'Hessian (symmetric matrix)',stH});
    map = colormap('gray'); colormap(flipud(map));
    colorbar; axis equal; axis tight; 
    
else
    figure; nr=3; nc=1;
    subplot(nr,nc,1);
    %pcolor(H); shading flat;
    imagesc(H);
    xlabel('Source index'); ylabel('Source index');
    title({'Hessian (symmetric matrix)',stH});
    map = colormap('gray'); colormap(flipud(map));
    colorbar; axis equal; axis tight; 

    subplot(nr,nc,2);
    %pcolor(Hstr); shading flat;
    imagesc(Hstr);
    xlabel('Source index'); ylabel('Source index');
    title('Hessian for structure (symmetric matrix)');
    map = colormap('gray'); colormap(flipud(map));
    colorbar; axis equal; axis tight; 

    subplot(nr,nc,3);
    %pcolor(Hsrc); shading flat;
    imagesc(Hsrc);
    xlabel('Source index'); ylabel('Source index');
    title('Hessian for source (diagonal matrix)');
    map = colormap('gray'); colormap(flipud(map));
    colorbar; axis equal; axis tight; 
    orient tall, wysiwyg
end

% spectrum of singular values
[U,S,V] = svd(H);
s = diag(S);
figure; semilogy(s,'.-','markersize',18);
xlabel('Index'); ylabel('Singular value');
title(sprintf('singular values of H, cond(H) = %.2e',cond(H)));

%=========================================================
