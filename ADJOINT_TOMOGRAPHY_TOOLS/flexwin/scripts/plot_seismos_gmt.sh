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
#t_end=200

# find number of windows
nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`

# default settings
paper="a4+"
orient="-P"
proj=X8/1.5
origin="-Y6"
shift="-Y-2.5"

# carl settings
#paper="letter"
#orient=" "
#proj=X10/1.8
#origin="-X0.5 -Y6"
#shift="-Y-2.5"

gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 14p LABEL_FONT_SIZE 16p
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
out=${basename}.seis.ps
outpdf=${basename}.seis.pdf

## for SCSN plots: 
#dkm=`grep DKM $info | awk '{print $NF}'`
#scsn_1=`echo $dkm | awk '{print -5+$1/7.5}'`
#scsn_2=`echo $dkm | awk '{print -15+$1/3.5}'`
#scsn_3=`echo $dkm | awk '{print 35+$1/3.1}'`

# set region
max=`grep PLOT_MAX ${obs} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

# do plot

psbasemap -R$region -J$proj -B${t_step}::.${label}:"Time (s)":/S -K $orient $origin > $out
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -m -W1 -G$paleblue -O -K >> $out
psxy $obs -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $syn -R$region -J$proj -H$nhead -W2/$red -O -K >> $out

#SCSN stuff
#psxy -R$region -J$proj -W2/${blue} -m -O -K << END >> $out
#$scsn_1 -$max
#$scsn_1 $max
#>
#$scsn_2 -$max
#$scsn_2 $max
#>
#$scsn_3 -$max
#$scsn_3 $max
#>
#END

#tail -$nwin $win | awk '{printf "> %f %f 8 0 0 LB 10p 0.5i l\n%s\n", $2,-m,substr($0,31,200)}' m=$max | pstext -R$region -J$proj -m -N -O -K -V >> $out
tail -$nwin $qual | awk '{printf "> %f %f 8 90 0 RT 10p 0.5i l\nCC=%.2f\ndT=%.2f\ndA=%.2f\n",$2,m,$5,$4,$6}' m=$max | pstext -R$region -J$proj -m -N -O -K >> $out
#tail -$nwin $qual | awk '{printf "> %f %f 8 90 0 RT 10p 0.5i l\nCC=%.2f\n", $2,m,$5,$4}' m=$max | pstext -R$region -J$proj -m -N -O -K >> $out
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 12 0 0 LT Seismograms
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# do plot for STA/LTA
#################################

nhead=8

# set region
max=`grep STA_LTA_MAX ${stalta} | awk '{print $NF}'`
region=$t_start/$t_end/0/$max

psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -m -W1 -G$paleblue -O -K >> $out

awk '{print $1, $2}' $stalta | psxy -R$region -J$proj -H$nhead -W2/$blue -O -K >> $out
awk '{print $1, $3}' $stalta | psxy -R$region -J$proj -H$nhead -W2ta/$blue -O -K >> $out

region=$t_start/$t_end/0/1
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> $t_start 0.9 12 0 0 LT 1 10 l 
@;$blue;STA/LTA @;;
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# plot envelopes with windows
#################################

nhead=4

# set region
max=`grep PLOT_MAX ${envobs} | awk '{print $NF}'`
region=$t_start/$t_end/0/$max

# do plot
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -m -W1 -G$paleblue -O -K >> $out
psxy $envobs -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $envsyn -R$region -J$proj -H$nhead -W2/$red -O -K >> $out
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 12 0 0 LT Envelopes
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -U/0/-0.75/$basename >> $out
#psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out

#bbox=`grep "%%BoundingBox: " $out | tail -1 `
#sed s/"$bbox"/"%"/ $out | sed s/"%%BoundingBox: (atend)"/"$bbox"/ > t1.$$
#epsffit -c -s 10 10 600 780 t1.$$ $out
#epstopdf $out
#rm t1.$$

ps2pdf $out $outpdf

#echo $out
#gv $out

#===================================================


