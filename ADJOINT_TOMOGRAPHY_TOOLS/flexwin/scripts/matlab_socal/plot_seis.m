
tmin = 10;
tmax = 100;

dir = '/net/denali/scratch1/carltape/svn/specfem/automeasure_data/socal_3Dsyn_dat/MEASURE/';

flabel = 'HEC.CI.BHZ';

[tdat,sdat] = textread([dir flabel '.obs'],'%f%f','headerlines',4);
[tsyn,ssyn] = textread([dir flabel '.syn'],'%f%f','headerlines',4);

tmin = min(tdat);
inds = find(and( tdat >= tmin, tdat <= tmax));

ymax = max(abs([sdat(inds) ; ssyn(inds)]));

figure; hold on;
plot(tdat,sdat,'b');
plot(tsyn,ssyn,'r--');
axis([tmin tmax -ymax ymax]); grid on;
legend('data','synthetics');
orient landscape, wysiwyg