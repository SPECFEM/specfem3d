% Demonstration of the trench-smoothing function "surf"
% in utils/dynamic_rupture/subductionzone_mesh/process_slab_rotate.py
% Its purpose is to modify a subduction fault geometry 
% by making the fault dip angle steeper at shallow depth, 
% to avoid elements with small angles at the trench.
%
% Change the two parameters below, to explore their effect on fault geometry
%
% J. P. Ampuero

DIP = 15; % dip angle (in degrees)
zcutTop = 8; % steepen the fault surface above this depth (in km)

%------------------------

% original planar dipping fault
x=[0.1:0.1:80];
z = -x*tand(DIP);

% shallow steepening (function "surf")
c=2/zcutTop;
znew = -log(exp(-c*z)-1)/c;

% fix the trench
xtrench = interp1(znew,x,0,'spline');
ii = find(x>xtrench, 1);
x(ii-1) = xtrench;
znew(ii-1) = 0;

% cut out everything above ground
znew = min( znew, 0 );

% figure
plot(-x,z, -x,znew);
xlabel('x (km)')
ylabel('z (km)')
axis equal
grid on
legend('original','modified')
