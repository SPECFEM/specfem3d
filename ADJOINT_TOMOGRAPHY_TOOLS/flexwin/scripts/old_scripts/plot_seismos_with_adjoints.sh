#!/bin/bash

basename=$1
out=${basename}.seis_adj.eps

obs=${basename}.obs
syn=${basename}.syn
envobs=${basename}.env.obs
envsyn=${basename}.env.syn
win=${basename}.win
qual=${basename}.win.qual
f1f2=${basename}.f1f2

seis=${basename}.seis.win

# find start and end time and a reasonable step 
t_start=`grep T_START ${obs} | awk '{print $NF}'`
t_end=`grep T_END ${obs} | awk '{print $NF}'`
t_step=`echo $t_end $t_start | awk '{print 50*int(($1-$2)/500)}'`

# find number of windows
nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`

if [ $nwin -gt 0 ] ; then

# find max min for adjoints
p_max=0
q_max=0

for file in ${seis}.*[0-9] ; do
  adj=`echo $file | sed s/seis.win/adj.win/`
  p_max_local=`grep P_MAX ${adj} | awk '{print $NF}'`
  q_max_local=`grep Q_MAX ${adj} | awk '{print $NF}'`
  
  p_max=`echo $p_max $p_max_local | awk '{if ( $2 > $1 ) print $2 ; else print $1}' `
  q_max=`echo $q_max $q_max_local | awk '{if ( $2 > $1 ) print $2 ; else print $1}' `
done

# if p_max or q_max are still = 0 (all adjoint sources are null,
# perfect fit), then set them to =1 for plotting
p_max=`echo $p_max | awk '{if ( $1 == 0.0 ) print 1.0 ; else print $1 }'`
q_max=`echo $q_max | awk '{if ( $1 == 0.0 ) print 1.0 ; else print $1 }'`

gmtset PAPER_MEDIA letter+ MEASURE_UNIT inch
proj=X8/3
nhead=4
red=255/0/0
paleblue=220/220/255

#################################
# plot seismograms with windows
#################################

# set output filename

# set region
max=`grep PLOT_MAX ${obs} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

# do plot
# SEISMOGRAM LINE
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -K -P -Y10 > $out
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -M -W1 -G$paleblue -O -K >> $out
psxy $obs -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $syn -R$region -J$proj -H$nhead -W2/$red -O -K >> $out
tail -$nwin $win | awk '{printf "> %f %f 8 0 0 LB 10p 0.5in l\n%s\n", $2,-m,substr($0,31,200)}' m=$max | pstext -R$region -J$proj -M -N -O -K >> $out
tail -$nwin $qual | awk '{printf "> %f %f 8 90 0 RT 10p 0.5in l\nF2=%.2f\ndT=%.2f\n", $2,m,$5,$6}' m=$max | pstext -R$region -J$proj -M -N -O -K >> $out
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 12 0 0 LT Seismograms
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out


# ADJOINT LINES
nhead=8
proj=X8/1.5

# travel time
region=$t_start/$t_end/-$p_max/$p_max
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K -Y-2.5 >> $out
psxy -R$region -J$proj -W2to -O -K << END >> $out
$t_start 0
$t_end 0
END
for file in ${seis}.*[0-9] ; do
  adj=`echo $file | sed s/seis.win/adj.win/`
  awk '{print $1, $2}' $adj | psxy -R$region -J$proj -H$nhead -W2 -O -K >> $out
done
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $p_max 12 0 0 LT dtau adjoint sources
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

# amplitude
region=$t_start/$t_end/-$q_max/$q_max
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K -Y-2.5 >> $out
psxy -R$region -J$proj -W2to -O -K << END >> $out
$t_start 0
$t_end 0
END
for file in ${seis}.*[0-9] ; do
  adj=`echo $file | sed s/seis.win/adj.win/`
  awk '{print $1, $3}' $adj | psxy -R$region -J$proj -H$nhead -W2 -O -K >> $out
done
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $q_max 12 0 0 LT dlnA adjoint sources
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out





# END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -U$basename >> $out

bbox=`grep "%%BoundingBox: " $out | tail -1 `
sed s/"$bbox"/"%"/ $out | sed s/"%%BoundingBox: (atend)"/"$bbox"/ > t1
epsffit -c -s 10 10 600 780 t1 $out

echo $out
#gv $out



fi
