#!/bin/csh -f

set depth = `grep "source depth" $1/MZZ/simulation.info |awk '{print $1}'`
set model = `grep Background $1/MZZ/mesh_params.h |awk '{print $5}'`
echo "taup depth: $depth, model: $model"

if ( $model == 'prem' || $model == 'iasp91' ) then
rm -f $2/taup_*_traveltime*
set num_rec = `wc -l MZZ/Data/receiver_pts.dat`
set epi_list = `tail -n $num_rec MZZ/Data/receiver_pts.dat | sed 's/999999/9/g' |  sed 's/000000/ /g' | awk '{print $1}'`
set depth_short = `echo $depth |sed 's/\./ /g' |awk '{print $1}'` 
echo $depth_short
set i = 0
foreach rec (${epi_list})
@ i++
echo "working on receiver" $rec
#taup_time -mod $model -h $depth -ph P,Pdiff -deg $rec 
#taup_time -mod $model -h $depth -ph S,P,p,s,PP,SS,ScS,PcP,PcS,PKIKP,PKIKS,Pdiff,Sdiff,PKP,PKS,pS,pP,PS,sS,SP,SSS,PPP,PSS,PSSS,P400P,P400S,P670P,P670S, -deg $rec | grep $depth_short   | awk '{print $4}' >! $2/taup_${model}_${depth}km_rec${i}_tt.dat

#taup_time -mod $model -h $depth -ph S,P,p,s,PP,SS,ScS,PcP,PcS,PKIKP,PKIKS,Pdiff,Sdiff,PKP,PKS,pS,pP,PS,sS,SP,SSS,PPP,PSS,PSSS,P400P,P400S,P670P,P670S, -deg $rec |grep $depth_short | awk '{print $8}' >! $2/taup_rec${i}_phase.dat
#cat  $2/taup_rec${i}_tt.dat |awk '{print $1 "  0.0"}' > $2/recfile_${i}_taup.dat
#paste  $2/taup_rec${i}_tt.dat $2/taup_rec${i}_phase.dat > $2/taup_rec${i}_tt_ph.dat
echo `taup_time -mod $model -h $depth -ph P -deg $rec | grep " P " | awk '{print $4}' |grep -v "==" `   $rec>> $2/taup_P_traveltime.dat
echo `taup_time -mod $model -h $depth -ph S -deg $rec | grep " S "  | awk '{print $4}' |grep -v "==" `  $rec >> $2/taup_S_traveltime.dat
echo `taup_time -mod $model -h $depth -ph PP -deg $rec | grep " PP " | awk '{print $4}' `  $rec>> $2/taup_PP_traveltime.dat
echo  `taup_time -mod $model -h $depth -ph PcP -deg $rec | grep " PcP "  | awk '{print $4}'` $rec >> $2/taup_PcP_traveltime.dat
echo  `taup_time -mod $model -h $depth -ph Pdiff -deg $rec | grep " Pdiff " | awk '{print $5}'` >> $2/taup_Pdiff_traveltime.dat
echo  `taup_time -mod $model -h $depth -ph Sdiff -deg $rec | grep " Sdiff " | awk '{print $5}'` >> $2/taup_Sdiff_traveltime.dat

end
endif

sort -n $2/taup_P_traveltime.dat |grep -v "==>" > $2/taup_P_traveltime2.dat
sort -n $2/taup_S_traveltime.dat |grep -v "==>" > $2/taup_S_traveltime2.dat
sort -n $2/taup_PP_traveltime.dat |grep -v "==>" > $2/taup_PP_traveltime2.dat
sort -n $2/taup_PcP_traveltime.dat |grep -v "==>" > $2/taup_PcP_traveltime2.dat
sort -n $2/taup_Pdiff_traveltime.dat |grep -v "==>" > $2/taup_Pdiff_traveltime2.dat
sort -n $2/taup_Sdiff_traveltime.dat |grep -v "==>" > $2/taup_Sdiff_traveltime2.dat

echo "Done with taup"
