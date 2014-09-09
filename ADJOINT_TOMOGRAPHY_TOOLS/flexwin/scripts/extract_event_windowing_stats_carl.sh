#!/bin/bash

basename=$1
figname=window_stats.fig

t0=t0
t1=t1
t2z=t2_Z
t2r=t2_R
t2t=t2_T
t3=t3

mincc=0.65
maxcc=1.01
ccstep=0.1

mindlnA=-1.5
maxdlnA=1.5
dlnAstep=0.5

minT=-10
maxT=10
Tstep=5

min_time=0
min_ddg=0

#CHT
rectxt=rectext

gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D LABEL_FONT_SIZE 12 ANOT_FONT_SIZE 10 HEADER_FONT_SIZE 12
red=255/0/0
black=0/0/0
grey=220/220/220
blue=0/0/255
green=00/200/00
paleblue=220/220/255

#################################
# extract window statistics
#################################

n_win_total=0
n_seis=0

max_ddg=0
max_time=0
for win in ${basename}/*.win.qual ; do
  echo $win
  nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`
  comp=`echo $win | awk -F"." '{print $(NF-2)}' `
  info=`echo $win | sed s/win.qual/info/`
  obs=`echo $win | sed s/win.qual/obs/`
  echo $comp $nwin >> $t1

  if [ $nwin -gt 0 ] ; then

    # extract single window information
    tail -$nwin $win | awk '{print c, $0}' c=$comp >> $t0

    # add to the number of seismograms and windows
    n_win_total=`echo $n_win_total $nwin | awk '{print $1+$2 }'`
    n_seis=`echo $n_seis | awk '{print $1+1 }'`

    # output the window recordsection information
    ddg=`grep DDG ${info} | awk '{print $NF}'`
    max_ddg=`echo $max_ddg $ddg | awk '{if ( $2 > $1) print $2 ; else print $1}'`
    t_end=`grep T_END ${obs} | awk '{print $NF}'`
    max_time=`echo $max_time $t_end | awk '{if ( $2 > $1) print $2 ; else print $1}'`

    c=`echo $comp | awk '{print substr($1,3,1)}'`
    if [ $c == "Z" ] ; then
      echo ">" >> $t2z
      echo "$min_time $ddg 0.0" >> $t2z
      tail -$nwin $win | awk '{printf "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", $2,d,0,$2,d,1,$3,d,1,$3,d,0}' d=$ddg >> $t2z
      echo "$t_end $ddg 0.0" >> $t2z
    elif [ $c == "R" ] ; then
      echo ">" >> $t2r
      echo "$min_time $ddg 0.0" >> $t2r
      tail -$nwin $win | awk '{printf "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", $2,d,0,$2,d,1,$3,d,1,$3,d,0}' d=$ddg >> $t2r
      echo "$t_end $ddg 0.0" >> $t2r
    elif [ $c == "T" ] ; then
      echo ">" >> $t2t
      echo "$min_time $ddg 0.0" >> $t2t
      tail -$nwin $win | awk '{printf "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", $2,d,0,$2,d,1,$3,d,1,$3,d,0}' d=$ddg >> $t2t
      echo "$t_end $ddg 0.0" >> $t2t
    fi

    # output the coverage information
    evla=`grep EVLA ${info} | awk '{printf "%.1f", $NF}'`
    evlo=`grep EVLO ${info} | awk '{printf "%.1f", $NF}'`
    stla=`grep STLA ${info} | awk '{printf "%.1f", $NF}'`
    stlo=`grep STLO ${info} | awk '{printf "%.1f", $NF}'`
    evdp=`grep EVDP ${info} | awk '{printf "%.1f", $NF}'`
    echo $evla $evlo $stla $stlo | awk '{printf ">\n%f %f\n%f %f\n",$1,$2,$3,$4}' >> $t3

  fi

done

timestep=`echo $min_time $max_time | awk '{print int(($2-$1)/6)}'`
degstep=`echo $min_ddg $max_ddg | awk '{print int(5*($2-$1)/25)}'`

#-------------
# CHT: make GMT text file for receivers
#iflag=0

max_time=165

#for win in ${basename}/*.win.qual ; do
#   strec=`basename $win | cut -d. -f1`
#   stnet=`basename $win | cut -d. -f2`
#   echo $strec $stnet | awk '{printf "%s %s\n",$1,$2}' >> $reclist
#done
#uniq_rec=`uniq reclist`

for win in ${basename}/*.win.qual ; do
  echo $win
  nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`
  info=`echo $win | sed s/win.qual/info/`

#  comp=`echo $win | awk -F"." '{print $(NF-3)}' `
#  suffix=`echo $win | awk -F"." '{print $(NF-2)}' `
  comp=`echo $win | awk -F"." '{print $(NF-2)}' `
  c=`echo $comp | awk '{print substr($1,3,1)}'`
  strec=`basename $win | cut -d. -f1`

  if [ $c == "R" ] ; then nwinR=$nwin ; fi
  if [ $c == "T" ] ; then nwinT=$nwin ; fi
  if [ $c == "Z" ] ; then nwinZ=$nwin ; fi

  # assume that there is always a Z component match
  if [ $c == "Z" ] ; then
    strec=`basename $win | cut -d. -f1`
    stnet=`basename $win | cut -d. -f2`
    ddg=`grep DDG ${info} | awk '{print $NF}'`      # arc distance
    #if [ `echo $iflag | awk '{print $1%2}'` == 0 ]; then
    #   tpos=$max_time
    #else 
    #   tpos=`echo $max_time $timestep | awk '{print $1+$2}'`
    #fi
    echo $max_time $ddg $strec $stnet $nwinZ $nwinR $nwinT | awk '{printf "%f %f %f %f %f %s %s %s (%s, %s, %s)\n",$1,$2,3,0,1,"LM",$3,$4,$5,$6,$7}' >> $rectxt
    #iflag=`echo $iflag | awk '{print $1+1}'`
    #echo $iflag

    #echo $strec $stnet | awk '{printf "%s %s\n",$1,$2}' >> $reclst
  fi

done
#-------------

#################################
# plot window statistics
#################################

# set output filename
out=${basename}/event_winstats.eps
out_ps=${basename}/event_winstats.ps
out_pdf=${basename}/event_winstats.pdf
short_basename=`echo $basename | awk -F"/" '{print $NF}'`

##evlo=`echo $evlo | awk '{printf "%.0f",  $1}'`
## set projection for map plotting
#proj=W$evlo/4
#echo $proj
#psbasemap -Rg -J$proj -Bg30/g30 -K -P -Y8 -V > $out
#pscoast -Rg -J$proj -Dl -A5000/0 -W1 -G$grey -O -K -V >> $out
#psxy $t3 -Rg -J$proj -: -m -W2 -O -K -V >> $out

# CHT
# set projection for map plotting
proj="-JM4"
bounds="-R-122/-114/32/37"
tick="-Ba1dWeSn"
echo $proj
psbasemap $bounds $proj $tick -K -Y7.5 -X0.8 -P -V > $out
pscoast $bounds $proj -Df -A100 -W1p -Na/1.0p -S150/200/255 -G200/255/150 -O -K -V >> $out
psxy $t3 $bounds $proj -: -m -W2 -O -K -V >> $out
psxy $proj $bounds -W1p -Sa0.2 -G250/0/0 -O -K >> $out <<EFX
$evlo $evla
EFX
psbasemap $bounds $proj $tick -O -K -V >> $out

proj=X3/2
pstext -R0/1/0/1 -J$proj -N  -O -K -X4.5 << END >> $out
0 1 14 0 0 LT $short_basename
0 0.80 12 0 0 LT Epicenter: ($evla, $evlo). 
0 0.65 12 0 0 LT Depth: $evdp km. 
0 0.50 12 0 0 LT Records with measurement windows: $n_seis. 
0 0.35 12 0 0 LT Measurement windows: $n_win_total.
END

proj=X1.75/1.5

## plot window number statistics
awk '{print $6}' $t0 | pshistogram  -J$proj -F -B:"CC":/WS -W0.02 -Lthin -G$grey -O -K -X-4.5 -Y-2.0 >> $out
awk '{print $5}' $t0 | pshistogram  -J$proj -F -B:"Tshift":/WS -W1 -Lthin -G$grey -O -K -X2.75 >> $out
awk '{print $7}' $t0 | pshistogram  -J$proj -F -B:"dlnA":/WS -W0.1 -Lthin -G$grey -O -X2.75 -K >> $out

# set projection for stats plotting
# plot Tshift vs CC 
region=$mincc/$maxcc/$minT/$maxT
psbasemap -R$region -J$proj -B${ccstep}:"CC":/${Tstep}:"Tshift / s":WS -O -K -Y-2.25 -X-5.5 >> $out
# plot the stats
#grep "Z" $t0 | awk '{print $6, $5}' | psxy -R$region -J$proj -Sl0.1/x -W2/$red -O -K -V >> $out
#grep "R" $t0 | awk '{print $6, $5}' | psxy -R$region -J$proj -Sl0.1/x -W2/$green -O -K -V >> $out
#grep "T" $t0 | awk '{print $6, $5}' | psxy -R$region -J$proj -Sl0.1/x -W2/$blue -O -K -V >> $out
awk '{print $6, $5}' $t0 | psxy -R$region -J$proj -Sl0.1/x -W2/$black -O -K -V >> $out
psbasemap -R$region -J$proj -B -O -K >> $out

# plot dlnA vs CC 
region=$mincc/$maxcc/$mindlnA/$maxdlnA
psbasemap -R$region -J$proj -B${ccstep}:"CC":/${dlnAstep}:"dlnA":WS -O -K -P -X2.75 >> $out
#grep "Z" $t0 | awk '{print $6, $7}' | psxy -R$region -J$proj -Sl0.1/x -W2/$red -O -K >> $out
#grep "R" $t0 | awk '{print $6, $7}' | psxy -R$region -J$proj -Sl0.1/x -W2/$green -O -K >> $out
#grep "T" $t0 | awk '{print $6, $7}' | psxy -R$region -J$proj -Sl0.1/x -W2/$blue -O -K >> $out
awk '{print $6, $7}' $t0 | psxy -R$region -J$proj -Sl0.1/x -W2/$black -O -K >> $out
psbasemap -R$region -J$proj -B -O -K >> $out

# plot dlnA vs Tshift 
region=$minT/$maxT/$mindlnA/$maxdlnA
psbasemap -R$region -J$proj -B${Tstep}:"Tshift":/${dlnAstep}:"dlnA":WS -O -K -P -X2.75 >> $out
#grep "Z" $t0 | awk '{print $5, $7}' | psxy -R$region -J$proj -Sl0.1/x -W2/$red -O -K >> $out
#grep "R" $t0 | awk '{print $5, $7}' | psxy -R$region -J$proj -Sl0.1/x -W2/$green -O -K >> $out
#grep "T" $t0 | awk '{print $5, $7}' | psxy -R$region -J$proj -Sl0.1/x -W2/$blue -O -K >> $out
awk '{print $5, $7}' $t0 | psxy -R$region -J$proj -Sl0.1/x -W2/$black -O -K >> $out
psbasemap -R$region -J$proj -B -O -K >> $out

## plot window number statistics
## NOTE: IF THERE IS A CONSTANT NUMBER OF WINDOWS FOR A COMPONENT, THEN PSHISTOGRAM CRASHES
grep Z $t1 | awk '{print $2}' | pshistogram  -J$proj -F -B:"Windows per Z trace":/WS -W1 -Lthin -G$red -O -K -X-5.5 -Y-2.25 >> $out
grep R $t1 | awk '{print $2}' | pshistogram  -J$proj -F -B:"Windows per R trace":/WS -W1 -Lthin -G$green -O -K -X2.75 >> $out
grep T $t1 | awk '{print $2}' | pshistogram  -J$proj -F -B:"Windows per T trace":/WS -W1 -Lthin -G$blue -O -X2.75 >> $out   # FINISH (no -K)

# CHT comment out
#echo "%!PS-Adobe-3.0" >t1
#echo "%%BoundingBox: 0 0 700 842" >> t1
#tail -n +3 $out >> t1
#mv t1 $out

echo $out

#gv $out

#-------------------------------

# plot window record-section

out2=${basename}/event_recordsection.eps
out2_ps=${basename}/event_recordsection.ps
out2_pdf=${basename}/event_recordsection.pdf
#gmtset PAPER_MEDIA letter MEASURE_UNIT inch

region=$min_time/$max_time/$min_ddg/$max_ddg
proj=X3.0/-6.8

psbasemap -R$region -J$proj -B${timestep}:"Time(s) -- Z":/${degstep}:"Distance (deg)":WN -K -X1.0 -Y0.75 -U/0/-0.25/${basename} > $out2
pswiggle $t2z -R$region -J$proj -Z20 -m -W1 -G$red -O -K >> $out2 

psbasemap -R$region -J$proj -B${timestep}:"Time(s) -- R":/${degstep}:"Distance (deg)":wN -O -K -X3.0 >> $out2
pswiggle $t2r -R$region -J$proj -Z20 -m -W1 -G$green -O -K >> $out2 

psbasemap -R$region -J$proj -B${timestep}:"Time(s) -- T":/${degstep}:"Distance (deg)":wN -O -K -X3.0 >> $out2
pswiggle $t2t -R$region -J$proj -Z20 -m -W1 -G$blue -O -K >> $out2 

pstext $rectxt -R$region -J$proj -N -O -G0 >> $out2   # FINISH (no -K)

echo $out2

#---------------------

#region=$min_time/$max_time/$min_ddg/$max_ddg
#proj=X5/-2

#psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WS -K -P -Y8 > $out2
#pswiggle $t2z -R$region -J$proj -Z20 -m -P -W1/$black -G$red -O -K >> $out2 
#psbasemap -R$region -J$proj -B:"Vertical component":/N -O -K >> $out2
#psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WS -O -K -Y-3.5 >> $out2
#pswiggle $t2r -R$region -J$proj -Z20 -m -P -W1/$black -G$green -O -K >> $out2 
#psbasemap -R$region -J$proj -B:"Radial component":/N -O -K >> $out2
#psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WS -O -K -Y-3.5 >> $out2
#pswiggle $t2t -R$region -J$proj -Z20 -m -P -W1/$black -G$blue -O -K >> $out2 
#psbasemap -R$region -J$proj -B:"Transverse component":/N -O >> $out2

#echo "%!PS-Adobe-3.0" >t1
#echo "%%BoundingBox: 0 0 500 842" >> t1
#tail -n +3 $out2 >> t1
#mv t1 $out2

#echo $out2
##gv $out

#fig_out_name=${basename}/event_composite.fig
#fig_out_eps=${basename}/event_composite.eps
#fig_out_pdf=${basename}/event_composite.pdf

#if [ -e $figname ] ; then
#  sed s!WINSTATS!$out! $figname | sed s!RECSECTION!$out2! > $fig_out_name
#  #fig2eps $fig_out_name
#  fig2pdf --nogv $fig_out_name
#  echo $fig_out_pdf
#fi

mv $rectxt $basename  # CHT: move pstext file to data directory

# make pdf files
cp $out ${out_ps}
ps2pdf ${out_ps} ${out_pdf}
cp $out2 ${out2_ps}
ps2pdf ${out2_ps} ${out2_pdf}

rm $t0 $t1 $t2z $t2r $t2t $t3

#===============================================
