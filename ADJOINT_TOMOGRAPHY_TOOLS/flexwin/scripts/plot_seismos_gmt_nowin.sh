#!/bin/bash

basename=$1
label=`basename $basename`

obs=${basename}.obs
syn=${basename}.syn
envobs=${basename}.env.obs
envsyn=${basename}.env.syn
win=${basename}.win
qual=${basename}.win.qual
stalta=${basename}.stalta
info=${basename}.info

# find start and end time and a reasonable step 
t_start=`grep T_START ${obs} | awk '{print $NF}'`
t_end=`grep T_END ${obs} | awk '{print $NF}'`
#t_step=`echo $t_end $t_start | awk '{print 50*int(($1-$2)/500)}'`
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/10)}'`

# CHT
#t_start=-20
#t_end=120
#t_end=200

# find number of windows
#nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`

gmtset PAPER_MEDIA a4+ MEASURE_UNIT inch HEADER_FONT_SIZE 14p LABEL_FONT_SIZE 16p
proj=X8/1.5
nhead=4
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255

#################################
# plot seismograms with windows
#################################

# set output filename
out=${basename}.seis_nowin.eps

# set region
max=`grep PLOT_MAX ${obs} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

# do plot
psbasemap -R$region -J$proj -B${t_step}::.${label}:"Time (s)":/S -K -P -Y8 > $out
psxy $syn -R$region -J$proj -H$nhead -W2 -O -K >> $out
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 12 0 0 LT Synthetic
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# plot envelopes with windows
#################################

nhead=4

# set region
max=`grep PLOT_MAX ${envobs} | awk '{print $NF}'`
region=$t_start/$t_end/0/$max
proj=X8/1.5

# do plot
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K -Y-2.5 >> $out
psxy $envsyn -R$region -J$proj -H$nhead -W2 -O -K >> $out
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 12 0 0 LT Envelope
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out


#################################
# do plot for STA/LTA
#################################

nhead=8

# set region
proj=X8/1.5
max=`grep STA_LTA_MAX ${stalta} | awk '{print $NF}'`
region=$t_start/$t_end/0/$max

psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K -Y-2.5 >> $out

awk '{print $1, $2}' $stalta | psxy -R$region -J$proj -H$nhead -W2 -O -K >> $out
awk '{print $1, $3}' $stalta | psxy -R$region -J$proj -H$nhead -W2ta -O -K >> $out

region=$t_start/$t_end/0/1
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start 0.9 12 0 0 LT STA:LTA
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out



bbox=`grep "%%BoundingBox: " $out | tail -1 `
sed s/"$bbox"/"%"/ $out | sed s/"%%BoundingBox: (atend)"/"$bbox"/ > t1.$$
epsffit -c -s 10 10 600 780 t1.$$ $out
epstopdf $out
rm t1.$$
#rm $out  # remove eps file

#echo $out
#gv $out

#===================================================


