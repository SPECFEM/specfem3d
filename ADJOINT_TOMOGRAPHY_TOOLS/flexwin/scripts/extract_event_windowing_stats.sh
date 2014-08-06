#!/bin/bash

if [ $# != 1 ]; then
  echo "extract_event_windowing_stats.sh basename(MEASURE)"; exit
fi

basename=$1
figname=window_stats.fig

t0=t0.$$
t1=t1.$$
t2z=t2_Z.$$
t2r=t2_R.$$
t2t=t2_T.$$
t3=t3.$$

mincc=0.65
maxcc=1.01
ccstep=0.1

mindlnA=-1.1
maxdlnA=1.1
dlnAstep=0.5

minT=-20
maxT=20
Tstep=5

min_time=0
min_ddg=0


gmtset PAPER_MEDIA a4+ MEASURE_UNIT inch
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


#################################
# plot window statistics
#################################

# set output filename
out=${basename}/event_winstats.eps
outpdf=${basename}/event_winstats.pdf
short_basename=`echo $basename | awk -F"/" '{print $NF}'`

#evlo=`echo $evlo | awk '{printf "%.0f",  $1}'`
# set projection for map plotting
proj=W$evlo/5
echo $proj
psbasemap -Rg -J$proj -Bg30/g30 -K -P -Y5 -V > $out
pscoast -Rg -J$proj -Dl -A5000/0 -W1 -G$grey -O -K -V >> $out
psxy $t3 -Rg -J$proj -: -m -W2 -O -K -V >> $out

proj=X3/2
pstext -R0/1/0/1 -J$proj -N  -O -K -X5.5 << END >> $out
0 1 18 0 0 LT $short_basename
0 0.80 12 0 0 LT Epicenter: ($evla, $evlo). 
0 0.65 12 0 0 LT Depth: $evdp km. 
0 0.50 12 0 0 LT N_rec : $n_seis. 
0 0.35 12 0 0 LT N_win : $n_win_total.
END

## plot window number statistics
proj=X2/2
#awk '{print $6}' $t0 | pshistogram  -J$proj -F -B:"CC":/WS -W0.02 -Lthin -G$grey -O -K -X-4.5 -Y-3.5 >> $out
awk '{print $6}' $t0 | pshistogram  -J$proj -F -B:"CC":/WS -W0.02 -Lthin -G$grey -O -K -X-5.5 -Y-3.5 >> $out
awk '{print $5}' $t0 | pshistogram  -J$proj -F -B:"@~t@~ / s":/WS -W1 -Lthin -G$grey -O -K -X3 >> $out
awk '{print $7}' $t0 | pshistogram  -J$proj -F -B:"@~D@~lnA":/WS -W0.1 -Lthin -G$grey -O -X3 >> $out


echo "%!PS-Adobe-3.0" >t1
echo "%%BoundingBox: 0 0 700 600" >> t1
tail -n +3 $out >> t1
mv t1 $out

echo $out
ps2pdf $out $outpdf
# plot window recordsection

out2=${basename}/event_recordsection.eps
out2pdf=${basename}/event_recordsection.pdf
gmtset PAPER_MEDIA a4+ MEASURE_UNIT inch

region=$min_time/$max_time/$min_ddg/$max_ddg
proj=X8/-2

psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WS -K -P -Y8 > $out2
pswiggle $t2z -R$region -J$proj -Z20 -m -P -W1/$black -G$red -O -K >> $out2 
psbasemap -R$region -J$proj -B:"Vertical component":/N -O -K >> $out2
psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WS -O -K -Y-3.5 >> $out2
pswiggle $t2r -R$region -J$proj -Z20 -m -P -W1/$black -G$green -O -K >> $out2 
psbasemap -R$region -J$proj -B:"Radial component":/N -O -K >> $out2
psbasemap -R$region -J$proj -B${timestep}:"Time(s)":/${degstep}:"Distance / degree":WS -O -K -Y-3.5 >> $out2
pswiggle $t2t -R$region -J$proj -Z20 -m -P -W1/$black -G$blue -O -K >> $out2 
psbasemap -R$region -J$proj -B:"Transverse component":/N -O >> $out2

echo "%!PS-Adobe-3.0" >t1
echo "%%BoundingBox: 0 0 700 842" >> t1
tail -n +3 $out2 >> t1
mv t1 $out2

echo $out2
ps2pdf $out2 $out2pdf

rm -f $t0 $t1 $t2z $t2r $t2t $t3

