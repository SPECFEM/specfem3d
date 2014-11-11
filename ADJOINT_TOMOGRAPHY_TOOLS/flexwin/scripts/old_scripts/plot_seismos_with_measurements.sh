#!/bin/bash

basename=$1
out=${basename}.seis_measure.eps

obs=${basename}.obs
syn=${basename}.syn
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
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -K -P -Y50 -X3 > $out
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -M -W1 -G$paleblue -O -K >> $out
psxy $obs -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $syn -R$region -J$proj -H$nhead -W2/$red -O -K >> $out
tail -$nwin $win | awk '{printf "> %f %f 8 0 0 LB 10p 0.5in l\n%s\n", $2,-m,substr($0,31,200)}' m=$max | pstext -R$region -J$proj -M -N -O -K >> $out
tail -$nwin $qual | awk '{printf "> %f %f 8 90 0 RT 10p 0.5in l\nF2=%.2f\ndT=%.2f\n", $2,m,$5,$6}' m=$max | pstext -R$region -J$proj -M -N -O -K >> $out
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 12 0 0 LT Seismograms
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out


echo Plotting measurements
# MEASUREMENTS
first_plot=1

for file in ${seis}.*[0-9] ; do

# window number
w_num=`echo $file | awk -F "." '{print $NF}'`
#echo $w_num First plot = $first_plot
# dtau dlnA files
dtau=`echo $file | sed s/seis.win/dtau/ `
dlnA=`echo $file | sed s/seis.win/dlnA/ `
adj=`echo $file | sed s/seis.win/adj.win/ `


# number of frequency points
npts_fr=`grep NPTS ${dtau} | awk '{print $NF}'`

if [ $npts_fr -gt 1 ] ; then
# find start and end frequency and a reasonable step 
f_start=`grep F_START ${dtau} | awk '{print $NF}'`
f_end=`grep F_END ${dtau} | awk '{print $NF}'`
f_step=`echo $f_end $f_start| awk '{print ($1-$2)/3}'`
f_step_log=`echo $f_step | awk '{print int(log($1)/log(10.0))-1}'`
f_step=`echo $f_step $f_step_log | awk '{print int($1/10^$2)*10^$2}'`

# find min and max value and a reasonable step 
v_min_dtau=`grep _MIN ${dtau} | awk '{print $NF}'`
v_max_dtau=`grep _MAX ${dtau} | awk '{print $NF}'`
v_step_dtau=`echo $v_max_dtau $v_min_dtau | awk '{print ($1-$2)/5}'`
log=`echo $v_step_dtau | awk '{print int(log($1)/log(10.0))-1}'`
v_step_dtau=`echo $v_step_dtau $log | awk '{print int($1/10^$2)*10^$2}'`

v_min_dlnA=`grep _MIN ${dlnA} | awk '{print $NF}'`
v_max_dlnA=`grep _MAX ${dlnA} | awk '{print $NF}'`
v_step_dlnA=`echo $v_max_dlnA $v_min_dlnA | awk '{print ($1-$2)/5}'`
log=`echo $v_step_dlnA | awk '{print int(log($1)/log(10.0))-1}'`
v_step_dlnA=`echo $v_step_dlnA $log | awk '{print int($1/10^$2)*10^$2}'`
fi

if [ $npts_fr -eq 1 ] ; then
# find values for single frequency measurement
f_xcorr=`grep F_START ${dtau} | awk '{print $NF}'`
v_dtau=`grep _MIN ${dtau} | awk '{printf "%.2e", $NF}'`
v_dlnA=`grep _MIN ${dlnA} | awk '{printf "%.2e", $NF}'`
fi

# find f1 and f2
f1=`grep F1 ${dtau} | awk '{printf "%.2f", $NF}'`
f2=`grep F2 ${dtau} | awk '{printf "%.2f", $NF}'`

# number of time points
npts_t=`grep NPTS ${file} | awk '{print $NF}'`

# find start and end time and a reasonable step 
t_start=`grep T_START ${file} | awk '{print $NF}'`
t_end=`grep T_END ${file} | awk '{print $NF}'`

t_step=`echo $t_end $t_start | awk '{print ($1-$2)/3}'`
log=`echo $t_step | awk '{print int(log($1)/log(10.0))-1}'`
t_step=`echo $t_step $log | awk '{print int($1/10^$2)*10^$2}'`

# find min and max value of seismogram
v_max_seis=`grep PLOT_MAX ${file} | awk '{print $NF}'`

# find min and max value of adjoint sources
p_max_adj=`grep P_MAX ${adj} | awk '{print $NF}'`
q_max_adj=`grep Q_MAX ${adj} | awk '{print $NF}'`

# adjoint sources *could* be zero everywhere (perfect fit)
# if so, then need to adapt plotting range
p_max_adj=`echo $p_max_adj | awk '{if ( $1 > 0 ) print $1; else print 1}'`
q_max_adj=`echo $q_max_adj | awk '{if ( $1 > 0 ) print $1; else print 1}'`

gmtset PAPER_MEDIA letter+ MEASURE_UNIT inch
proj=X1.5/1.5
nhead_freq=7
nhead_seis=8
nhead_adj=8
red=255/0/0
black=0/0/0

#################################
# plot results
#################################

gmtset BASEMAP_TYPE plain PAPER_MEDIA letter+ LABEL_FONT_SIZE 10 ANNOT_FONT_SIZE 10 HEADER_FONT_SIZE 14

if [ $npts_fr -gt 0 ] ; then


# set region for seismograms
region=$t_start/$t_end/-$v_max_seis/$v_max_seis

# Plot synthetic vs observed before measurement
if [ $first_plot -eq 1 ] ; then
  psbasemap -R$region -J$proj -B${t_step}:"Time / s"::."Before":/S -O -K -Y-2.8  >> $out
else
  psbasemap -R$region -J$proj -B${t_step}:"Time / s":/S -O -K -Y-2.5 -X-6  >> $out
fi
psxy -R$region -J$proj -W1ta -O -K << END >> $out
$t_start 0
$t_end 0
END
awk '{print $1, $3}' $file | psxy -R$region -J$proj -H$nhead_seis -W2/$black -O -K >> $out
awk '{print $1, $2}' $file | psxy -R$region -J$proj -H$nhead_seis -W2/$red -O -K >> $out
pstext -R$region -J$proj -D0/-0.1 -N -O -K  <<  END >> $out
$t_start $v_max_seis 10 0 1 LT ($w_num)
END
psbasemap -R$region -J$proj -BS -O -K >>  $out

# Plot reconstructed observed vs observed after measurement
if [ $first_plot -eq 1 ] ; then
psbasemap -R$region -J$proj -B${t_step}:"Time / s"::."After":/S -O -K -P -X2 >> $out
else
psbasemap -R$region -J$proj -B${t_step}:"Time / s":/S -O -K -P -X2 >> $out
fi
psxy -R$region -J$proj -W1ta -O -K << END >> $out
$t_start 0
$t_end 0
END
awk '{print $1, $3}' $file | psxy -R$region -J$proj -H$nhead_seis -W2/$black -O -K >> $out
awk '{print $1, $4}' $file | psxy -R$region -J$proj -H$nhead_seis -W2/$red -O -K >> $out
pstext -R$region -J$proj -D-0.20/0 -N -O -K <<  END >> $out
$t_start 0 10 90 0 CB (F1=$f1, F2=$f2)
END
psbasemap -R$region -J$proj -BS -O -K >>  $out


###################################
# PLOT FREQUENCY DEPENDENT MEASUREMENT
###################################

if [ $npts_fr -gt 1 ] ; then

# reset region for dtau
# turn all frequencies into mHz
f_start=`echo $f_start | awk '{print $1*1000}'`
f_end=`echo $f_end | awk '{print $1*1000}'`
f_step=`echo $f_step | awk '{print $1*1000}'`
f_step_log=`echo $f_step | awk '{print int(log($1)/log(10.0))-1}'`
f_step=`echo $f_step $f_step_log | awk '{print int($1/10^$2)*10^$2}'`
region=$f_start/$f_end/$v_min_dtau/$v_max_dtau

# do plot
if [ $first_plot -eq 1 ] ; then
psbasemap -R$region -J$proj -B${f_step}:."@~dt@~ / s"::"Frequency / mHz":/${v_step_dtau}SW -O -K  -X2 >> $out
else
psbasemap -R$region -J$proj -B${f_step}:"Frequency / mHz":/${v_step_dtau}SW -O -K  -X2 >> $out
fi
awk '{print $1*1000, $2, $3}' $dtau | psxy -R$region -J$proj -H$nhead_freq -Ey -W2 -O -K >> $out
awk '{print $1*1000, $2}' $dtau | psxy -R$region -J$proj -H$nhead_freq -W2 -O -K >> $out
psbasemap -R$region -J$proj -BSW -O -K >> $out

# reset region for dlnA
region=$f_start/$f_end/$v_min_dlnA/$v_max_dlnA

# do plot
if [ $first_plot -eq 1 ] ; then
psbasemap -R$region -J$proj -B${f_step}:"Frequency / mHz":/:."@~d@~lnA":${v_step_dlnA}SW -O -K -X2 >> $out
else
psbasemap -R$region -J$proj -B${f_step}:"Frequency / mHz":/${v_step_dlnA}SW -O -K -X2 >> $out
fi
awk '{print $1*1000, $2, $3}' $dlnA | psxy -R$region -J$proj -H$nhead_freq -Ey -W2 -O -K >> $out
awk '{print $1*1000, $2}' $dlnA | psxy -R$region -J$proj -H$nhead_freq -W2 -O -K >> $out
psbasemap -R$region -J$proj -BSW -O -K >> $out

fi

###################################
# END FREQUENCY DEPENDENT MEASUREMENT
###################################

###################################
# PLOT XCORR MEASUREMENT
###################################
if [ $npts_fr -eq 1 ] ; then
region=0/10/-1/1
pstext -R$region -J$proj -O -K -X2 << END >> $out
5 0 12 0 0 CM @~dt@~=${v_dtau}s
END
pstext -R$region -J$proj -O -K -X2 << END >> $out
5 0 12 0 0 CM @~d@~lnA=${v_dlnA}
END
fi


# if we did a plot, and it was a first plot, then set first_plot flag to 0
if [ $first_plot -eq 1 ] ; then 
first_plot=0 
fi


fi # end if n_freq gt 0

done



# END
region=0/1/-1/1
#pstext -R$region -J$proj -O -X-6 -Y-0.2 -U$basename << END >> $out
pstext -R$region -J$proj -O -X-6 -Y-0.2  << END >> $out
END


bbox=`grep "%%BoundingBox: " $out | tail -1 `
sed s/"$bbox"/"%"/ $out | sed s/"%%BoundingBox: (atend)"/"$bbox"/ > t1
epsffit -c -s 10 10 600 780 t1 $out

echo $out
#gv $out

fi


