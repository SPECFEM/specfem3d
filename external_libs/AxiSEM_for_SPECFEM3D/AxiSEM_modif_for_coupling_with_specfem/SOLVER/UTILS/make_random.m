function make_random(rmin, rmax, corr_length, maxval)

dx = corr_length / 2;
[x, y] = meshgrid((1:dx:rmax), (-rmax:dx:rmax));
r = sqrt(x.^2 + y.^2);

rand_ac = exp(-(x.^2 + y.^2)/corr_length^2);


rand_ft = fft2(rand_ac);
rand_field_ft = rand(size(rand_ft))*2 - 1;
% figure;
% imagesc(rand_ac);
% figure
rand_field = real(ifft2(rand_ft .* exp(2*pi*1i*rand_field_ft)));
% imagesc(rand_field)
% imagesc(r)

rand_field((r<rmin) | (r>rmax)) = 0.0;
rand_field = rand_field / max(max(rand_field)) * maxval;
figure;
imagesc(rand_field)


theta = atan2(y,x) * 180 / pi + 90;

% figure;
% imagesc(theta)

fid = fopen('random.het', 'w');
fprintf(fid, '%d\n', length(find((r>rmin-10) & (r<rmax+10))));
for ix = 1:size(rand_field, 1)
    for iy = 1:size(rand_field, 2)
        if (r(ix, iy) > rmin-10 && r(ix, iy) < rmax+10)
            fprintf(fid, '%f %f %f %f %f\n', r(ix, iy), theta(ix, iy), ...
                    rand_field(ix, iy), rand_field(ix, iy), rand_field(ix, iy));
        end 
    end 
end
        
fclose(fid);

end 