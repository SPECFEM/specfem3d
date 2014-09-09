#!/bin/bash

basename=$1

t0=extract_event_t0
t1=extract_event_t1
t2z=extract_event_t2_Z
t2r=extract_event_t2_R
t2t=extract_event_t2_T
t3=extract_event_t3
rectxt=rectext            # CHT
#reclst=reclist            # CHT

minf1=0.49
maxf1=1.01
minf2=0.49
maxf2=1.01
f2step=0.1

mindlnA=-1.1
maxdlnA=1.1
dlnAstep=0.2

minT=-10.5
maxT=10.5
Tstep=2

min_time=0
min_ddg=0

#timestep=1000
#degstep=20


gmtset PAPER_MEDIA letter+ MEASURE_UNIT inch ANNOT_FONT_SIZE_SECONDARY 8 LABEL_FONT_SIZE 10 HEADER_FONT_SIZE 10 ANNOT_FONT_SIZE 10
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255

#################################
# extract window statistics
#################################

if [ -e $t0 ] ; then ( rm $t0 ) ; fi
if [ -e $t1 ] ; then ( rm $t1 ) ; fi
if [ -e $t2z ] ; then ( rm $t2z ) ; fi
if [ -e $t2r ] ; then ( rm $t2r ) ; fi
if [ -e $t2t ] ; then ( rm $t2t ) ; fi
if [ -e $t3 ] ; then ( rm $t3 ) ; fi
if [ -e $rectxt ] ; then ( rm $rectxt ) ; fi    # CHT

n_win_total=0
n_seis=0
n_rej_total=0
n_rej_f2=0
n_rej_dlnA=0
n_rej_Tshift=0


max_ddg=0
max_time=0
for win in ${basename}/*.win.qual ; do
  echo $win
  nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`
#  comp=`echo $win | awk -F"." '{print $(NF-3)}' `
  comp=`echo $win | awk -F"." '{print $(NF-2)}' `
  info=`echo $win | sed s/win.qual/info/`
  obs=`echo $win | sed s/win.qual/obs/`
  echo $comp $nwin >> $t1

  if [ $nwin -gt 0 ] ; then

  rej_total=`grep NUM_REJ_TOTAL ${info} | awk '{print $NF}'`
  rej_f2=`grep NUM_REJ_F2 ${info} | awk '{print $NF}'`
  rej_dlnA=`grep NUM_REJ_DLNA ${info} | awk '{print $NF}'`
  rej_Tshift=`grep NUM_REJ_TSHIFT ${info} | awk '{print $NF}'`
  
  n_rej_total=`echo $n_rej_total $rej_total | awk '{print $1+$2 }'`
  n_rej_f2=`echo $n_rej_f2 $rej_f2 | awk '{print $1+$2 }'`
  n_rej_dlnA=`echo $n_rej_dlnA $rej_dlnA | awk '{print $1+$2 }'`
  n_rej_Tshift=`echo $n_rej_Tshift $rej_Tshift | awk '{print $1+$2 }'`

    # extract single window information
    tail -$nwin $win | awk '{print c, $0}' c=$comp >> $t0

    # add to the number of seismograms and windows
    n_win_total=`echo $n_win_total $nwin | awk '{print $1+$2 }'`
    n_seis=`echo $n_seis | awk '{print $1+1 }'`

    # output the window record section information
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

timestep=`echo $min_time $max_time | awk '{print int(($2-$1)/4)}'`
#degstep=`echo $min_ddg $max_ddg | awk '{print int(5*($2-$1)/50)}'`
degstep=`echo $min_ddg $max_ddg | awk '{print int(10*($2-$1))/(10*10)}'`

max_ddg=`echo $max_ddg $degstep | awk '{print $1+$2/2}'`

echo $min_ddg $max_ddg $degstep

#-------------
# CHT: make GMT text file for receivers
#iflag=0
for win in ${basename}/*.win.qual ; do
  echo $win
  nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`
  info=`echo $win | sed s/win.qual/info/`

#  comp=`echo $win | awk -F"." '{print $(NF-3)}' `
#  suffix=`echo $win | awk -F"." '{print $(NF-2)}' `
  comp=`echo $win | awk -F"." '{print $(NF-2)}' `
  c=`echo $comp | awk '{print substr($1,3,1)}'`

  if [ $c == "R" ] ; then nwinR=$nwin ; fi
  if [ $c == "T" ] ; then nwinT=$nwin ; fi
  if [ $c == "Z" ] ; then
    nwinZ=$nwin
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

#echo "we are done"
#exit

#################################
# plot window statistics
#################################

# set output filename
out=${basename}/event_winstats.eps

evlo=`echo $evlo | awk '{printf "%.0f",  $1}'`
# set projectio for map plotting
proj="-JM1.5"
bounds="-R108/152/18/51"
tick="-Ba10dWeSn"
echo $proj
psbasemap $bounds $proj $tick -K -Y12 -X0.8 -P -V > $out
pscoast $bounds $proj -Dh -A100 -W0.5 -S150/200/255 -W2 -G200/255/150  -O -K -V >> $out
psxy $t3 $bounds $proj -: -m -W2 -O -K -V >> $out
psxy $proj $bounds -W1 -Sa0.2 -G250/0/0 -O -K -: >> $out <<EFX
$evla $evlo
EFX

proj=X3/2
pstext -R0/1/0/1 -J$proj -N  -O -K -X2 -Y-0.6<< END >> $out
0 1 14 0 0 LT $basename
0 0.80 12 0 0 LT Epicenter: ($evla, $evlo). 
0 0.65 12 0 0 LT Depth: $evdp km. 
0 0.50 12 0 0 LT Records with measurement windows: $n_seis. 
0 0.35 12 0 0 LT Measurement windows: $n_win_total.
#0 0.20 12 0 0 LT Rejection (Total: F2, dlnA, Tshift) : $n_rej_total : $n_rej_f2, $n_rej_dlnA, $n_rej_Tshift.
END

# NOTE: HISTOGRAMS CAN CAUSE PROBLEMS IF THERE ARE TOO FEW POINTS

# plot window number statistics
proj=X1.4/1.4
grep Z $t1 | awk '{print $2}' | pshistogram  -J$proj -C -B:"Windows per Z trace":/WS -W1 -Lthin -G$red -O -K -X-2.1 -Y-1.4 >> $out
grep R $t1 | awk '{print $2}' | pshistogram  -J$proj -C -B:"Windows per R trace":/WS -W1 -Lthin -G$green -O -K -X1.9 >> $out
grep T $t1 | awk '{print $2}' | pshistogram  -J$proj -C -B:"Windows per T trace":/WS -W1 -Lthin -G$blue -O -X1.9 -K >> $out

# plot window record section

region=$min_time/$max_time/$min_ddg/$max_ddg
proj=X1.4/-2.5

psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WN -O -K -X-3.8 -Y-3.7 -U/0/-0.25/${basename} >> $out
pswiggle $t2z -R$region -J$proj -Z20 -m -W1 -G$red -O -K >> $out 
psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}WN -O -K -X1.9>> $out
pswiggle $t2r -R$region -J$proj -Z20 -m -W1 -G$green -O -K >> $out 
psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}WN -O -K -X1.9 >> $out
pswiggle $t2t -R$region -J$proj -Z20 -m -W1 -G$blue -O >> $out 

bbox=`grep "%%BoundingBox: " $out | tail -1 `
sed s/"$bbox"/"%"/ $out | sed s/"%%BoundingBox: (atend)"/"$bbox"/ > t1
epsffit -c -s 10 10 600 800 t1 $out

echo $out
#ggv $out &

#--------------------
# CHT: repeat the window record section in another figure
# Landscape: -X(to the right) -Y(up)

# set output filename
out=${basename}/event_winstats_rs.ps

proj=X3.0/-6.8    # negative sign to flip the y-axis

psbasemap -R$region -J$proj -B${timestep}:"Time(s) -- Z":/${degstep}:"Distance (deg)":WN -K -X1.0 -Y0.75 -U/0/-0.25/${basename} > $out
pswiggle $t2z -R$region -J$proj -Z20 -m -W1 -G$red -O -K >> $out 

psbasemap -R$region -J$proj -B${timestep}:"Time(s) -- R":/${degstep}:"Distance (deg)":wN -O -K -X3.0 >> $out
pswiggle $t2r -R$region -J$proj -Z20 -m -W1 -G$green -O -K >> $out 

psbasemap -R$region -J$proj -B${timestep}:"Time(s) -- T":/${degstep}:"Distance (deg)":wN -O -K -X3.0 >> $out
pswiggle $t2t -R$region -J$proj -Z20 -m -W1 -G$blue -O -K >> $out 

pstext $rectxt -R$region -J$proj -N -O -G0 >> $out   # FINISH (no -K)

#bbox=`grep "%%BoundingBox: " $out | tail -1 `
#sed s/"$bbox"/"%"/ $out | sed s/"%%BoundingBox: (atend)"/"$bbox"/ > t1
#epsffit -c -s 10 10 600 800 t1 $out

echo $out
#ggv $out

#--------------------

rm $t0 $t1 $t2z $t2r $t2t $t3

mv $rectxt $basename  # move pstext file to data directory
#mv $reclst $basename  # move station list to data directory
