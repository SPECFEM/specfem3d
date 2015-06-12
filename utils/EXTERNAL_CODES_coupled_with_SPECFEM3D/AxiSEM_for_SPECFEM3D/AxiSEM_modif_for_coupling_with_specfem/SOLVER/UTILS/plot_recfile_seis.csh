#!/bin/csh -f

rm -f disp_*.ps disp_*.gif 

set list = `ls Data/*_disp.dat`
cat Data/seisepicenter1.dat |awk '{print $1}' > timeseries.dat
set mint = `head -n 1 timeseries.dat`; set maxt = `tail -n 1 timeseries.dat`
set comp1 = 's'; set comp2 = 'ph'; set comp3 = 'z'
set rot_rec = `grep rot_rec ../inparam |awk '{print $1}'`

if ( $rot_rec == 'T') then 
    set comp1 = 'r'; set comp2 = 'ph'; set comp3 = 'th'
else 
    set comp1 = 's'; set comp2 = 'ph'; set comp3 = 'z'
endif

echo "displacement components:" $comp1,$comp2,$comp3

set recname = `cat Data/receiver_names.dat`

set ind = 0

foreach file ( $list )

@ ind++

set i = `echo $file |sed 's/\.dat/ /' |sed 's/\./ /' |awk '{print $2}'`

set outfile1 = $recname[$ind]"_disp_"$comp1
set outfile3 = $recname[$ind]"_disp_"$comp3
if ( -f Data/seisepicenter2.dat) then 
    set outfile2 = $recname[$ind]"_disp_"$comp2
endif

set outfile1 = `echo $file"_"$comp1`
set outfile2 = `echo $file"_"$comp2`
set outfile3 = `echo $file"_"$comp3`

#set filen = `echo '"'"$file"'"'`

cat $file |awk '{print $1}' >poo 
awk '{str=$1; getline < "poo"; print str,$1}' timeseries.dat > $outfile1.dat; rm -f poo

if ( -f Data/seisepicenter2.dat) then 
#    echo "multipole source including phi component"
    echo 'working on' $outfile1,$outfile2,$outfile3
    cat $file |awk '{print $2}' > poo 
    awk '{str=$1; getline < "poo"; print str,$1}' timeseries.dat > $outfile2.dat; rm -f poo
    cat $file |awk '{print $3}' > poo 
    awk '{str=$1; getline < "poo"; print str,$1}' timeseries.dat > $outfile3.dat; rm -f poo

else # monopole, no phi component, second column is z/theta component
    cat $file |awk '{print $2}' > poo 
    awk '{str=$1; getline < "poo"; print str,$1}' timeseries.dat > $outfile3.dat; rm -f poo
	echo 'working on' $outfile1,$outfile3
endif

set minu1 = `minmax -C $outfile1.dat |awk '{print $3}'`
set maxu1 = `minmax -C $outfile1.dat  |awk '{print $4}'`
set minu3 = `minmax -C $outfile3.dat |awk '{print $3}'`
set maxu3 = `minmax -C $outfile3.dat |awk '{print $4}'`
if ( -f Data/seisepicenter2.dat) then
    set minu2 = `minmax -C $outfile2.dat |awk '{print $3}'`
    set maxu2 = `minmax -C $outfile2.dat |awk '{print $4}'`
endif

xmgrace -hdevice PostScript -hardcopy $outfile1.dat -legend load -world $mint $minu1 $maxt $maxu1 -printfile $outfile1.ps
convert -rotate 90 $outfile1.ps $outfile1.gif; 

xmgrace -hdevice PostScript -hardcopy $outfile3.dat -legend load -world $mint $minu3 $maxt $maxu3 -printfile $outfile3.ps
convert -rotate 90 $outfile3.ps $outfile3.gif; 

if ( -f Data/seisepicenter2.dat) then
    xmgrace -hdevice PostScript -hardcopy $outfile2.dat -legend load -world $mint $minu2 $maxt $maxu2 -printfile $outfile2.ps
    convert -rotate 90 $outfile2.ps $outfile2.gif; 
endif
gzip $file

gzip $outfile1.dat $outfile3.dat
if ( -f Data/seisepicenter2.dat) then
gzip $outfile2.dat
endif

end



