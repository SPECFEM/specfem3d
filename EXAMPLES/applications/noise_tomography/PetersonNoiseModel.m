function [accel, veloc, displ] = PetersonNoiseModel(Period,s)

if s=='NLNM'       % Peterson's New Low Noise Model
    P=[  0.10   0.17   0.40   0.80   1.24   2.40     4.30 ...
         5.00   6.00  10.00  12.00  15.60  21.90    31.60 ...
        45.00  70.00 101.00 154.00 328.00 600.00 10000.00 ];
    A=[ -162.36 -166.70 -170.00 -166.40 -168.60 -159.98 -141.10 ...
         -71.36  -97.26 -132.18 -205.27  -37.65 -114.37 -160.58 ...
        -187.50 -216.47 -185.00 -168.34 -217.43 -258.28 -346.88 ];
    B=[   5.64   0.00  -8.30  28.90   52.48  29.81   0.00 ...
        -99.77 -66.49 -31.57  36.16 -104.33 -47.10 -16.28 ...
          0.00  15.70   0.00  -7.61   11.90  26.60  48.75 ];
elseif s=='NHNM'   % Peterson's New High Noise Model
    P=[ 0.10 0.22  0.32  0.80   3.80 4.60 ...
        6.30 7.90 15.40 20.00   354.80 ];
    A=[ -108.73 -150.34 -122.31 -116.85 -108.48 -74.66 ...
           0.66  -93.37   73.54 -151.52 -206.66 ];
    B=[  -17.23 -80.50  -23.87 32.51 18.08 -32.95 ...
        -127.18 -22.42 -162.98 10.01 31.63 ];
end

N=length(P);

l=1;
% Period is less than P(1) (minimum period defined in the model)
if Period < P(1)
    accel=A(1)+B(1)*log10(P(1));
    veloc=accel+20.0*log10(P(1)/2/pi);
    displ=accel+20.0*log10(P(1)^2/4/pi^2);
end
% Period is between P(1) and P(N)
while l<=N-1
    if Period >= P(l) && Period < P(l+1)
        accel=A(l)+B(l)*log10(Period);
        veloc=accel+20.0*log10(Period/2/pi);
        displ=accel+20.0*log10(Period^2/4/pi^2);
        break;
    else
        l=l+1;
    end
end
% Period is larger than P(N) and less than 1e5 (maximum period defined in the model)
if Period>=P(N) && Period < 1e5
    accel=A(N)+B(N)*log10(Period);
    veloc=accel+20.0*log10(Period/2/pi);
    displ=accel+20.0*log10(Period^2/4/pi^2);
elseif Period > 1e5
    accel=A(N)+B(N)*log10(1e5);
    veloc=accel+20.0*log10(1e5/2/pi);
    displ=accel+20.0*log10(1e5^2/4/pi^2);
end
