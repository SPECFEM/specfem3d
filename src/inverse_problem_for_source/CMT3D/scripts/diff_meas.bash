# diff measurement files
for freq in 2 3 6; do
for evid in `cat event_list`; do
 res1=`head -n 1 $evid/MEASUREMENT_WINDOWS_${evid}_T00${freq}_T030_m12 `
 res2=`head -n 1 event_info/MEASUREMENT_WINDOWS_${evid}_T00${freq}_T030_m12`
 if [ $res1 -ne $res2 ]; then
 echo $evid -- $freq --- $res1 --- $res2
 fi
 done
done

