#!/bin/sh

# creates a simple, formatted tomography model with a constant velocity gradient 
# for a block model with dimensions 134000 x 134000 x 60000

# origin points
ORIG_X=0.
ORIG_Y=0. 
ORIG_Z=0. 

# end points
END_X=134000. 
END_Y=134000. 
END_Z=-60000.  # depth in negative z-direction

# spacing of given tomography points
SPACING_X=2000. 
SPACING_Y=2000. 
SPACING_Z=-2000.

# number of cell increments
NX=68
NY=68
NZ=31

# min/max values
VP_MIN=2500. 
VP_MAX=8500. 
VS_MIN=1500.
VS_MAX=7500. 
RHO_MIN=1500. 
RHO_MAX=1500.


# header info
echo "creating header info..."

echo "$ORIG_X $ORIG_Y $ORIG_Z $END_X $END_Y $END_Z  " > tmp.xyz
echo "$SPACING_X $SPACING_Y $SPACING_Z " >> tmp.xyz
echo "$NX $NY $NZ " >> tmp.xyz
echo "$VP_MIN $VP_MAX $VS_MIN $VS_MAX $RHO_MIN $RHO_MAX" >> tmp.xyz

# velocity gradient
GRADIENT=0.1

# adds point location and velocity model values
echo "adding model values..."

# format: lists first all x, then y, then z
echo "1" | awk '{ for(k=0;k<NZ;k++){ for(j=0;j<NY;j++){for(i=0;i<NX;i++){ x=i*SPACING_X;y=j*SPACING_Y;z=k*SPACING_Z;vp=VP_MIN+GRADIENT*(-z);vs=VS_MIN + GRADIENT*(-z); rho=RHO_MIN;print x,y,z,vp,vs,rho }}} }' \
            NX=$NX NY=$NY NZ=$NZ SPACING_X=$SPACING_X SPACING_Y=$SPACING_Y SPACING_Z=$SPACING_Z VP_MIN=$VP_MIN VS_MIN=$VS_MIN RHO_MIN=$RHO_MIN GRADIENT=$GRADIENT >> tmp.xyz

mv tmp.xyz tomography_model.xyz

echo "created file: tomography_model.xyz"
