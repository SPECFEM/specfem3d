#!/bin/bash
#
# creates adjoint traces
#
###############################################
# USER PARAMETERS

# station
#station="RS52S"
station="GRO1S"
network="EX"

# simulation setup
# minimum period
Tmin=0.3255

# non-zero trace component (1=X/2=Y/3=Z)
comp=2

band=HX

################################################

sta="$network.$station"

# time steps
NSTEP=`grep ^NSTEP ../DATA/Par_file | cut -d = -f 2 | sed "s/^[ \t]*//"`
dt=`grep ^DT ../DATA/Par_file | cut -d = -f 2 | sed "s/^[ \t]*//"`

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "creating adjoint source"
echo "  adjoint station: $sta"
echo "  component: $comp (1=X/2=Y/3=Z)"
echo "  length: $NSTEP"
echo "  time step size: $dt"
echo "  minimum period: $Tmin"
echo

# adjoint trace
# ricker wavelet
PI=3.141592653589793
f=`echo "1.0/$Tmin" | bc -l`
supp=`echo "($Tmin/$dt + 0.5)/1" | bc`

echo "ricker frequency: $f"
echo "ricker time step support: $supp"

tmp=trace.adj

echo "1" | awk '{a=PI**2 * f**2; \
                 for(i=0;i<NSTEP;i++){ \
                   t=(NSTEP-i)*dt - supp*dt; \
                   if (-a*t**2 < -500.0){ val=0.0;}else{ val=(1.0-2*a*t**2)*exp(-a*t**2);} \
                   print i*dt,"\t",val,"\t",t; \
                  }\
                 }' NSTEP=$NSTEP dt=$dt PI=$PI f=$f supp=$supp > $tmp
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# sets components
if [ "$comp" == "1" ]; then
  awk '{print $1,$2}'  $tmp > $sta.${band}X.adj
  awk '{print $1,"0.0"}' $tmp > $sta.${band}Y.adj
  awk '{print $1,"0.0"}' $tmp > $sta.${band}Z.adj
elif [ "$comp" == "2" ]; then
  awk '{print $1,"0.0"}' $tmp > $sta.${band}X.adj
  awk '{print $1,$2}'  $tmp > $sta.${band}Y.adj
  awk '{print $1,"0.0"}' $tmp > $sta.${band}Z.adj
elif [ "$comp" == "3" ]; then
  awk '{print $1,"0.0"}' $tmp > $sta.${band}X.adj
  awk '{print $1,"0.0"}' $tmp > $sta.${band}Y.adj
  awk '{print $1,$2}'  $tmp > $sta.${band}Z.adj
else
  echo "wrong component $comp"
  exit 1
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"


