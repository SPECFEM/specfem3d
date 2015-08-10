#!/bin/csh -f

echo ""
echo ""
echo ">>>>>>>>> AXISEM POST PROCESSING <<<<<<<<<"
echo ""
echo ""

if ( ${#argv} < 4 ) then
    if ( $1 == '-h' ) then
        echo " Input for post-processing is retrieved from input file:"
        echo "   1) param_post_processing : contains everything related to summation, convolution, receiver components"
        echo "   2) param_snaps : if wavefield snaps were computed, this file defines the 3D geometry"
        echo "   3) <simdir>/simulation.info: parameters from the simulation"
        echo
        exit
    endif
endif

set homedir = $PWD
set simtype = `grep "SIMULATION_TYPE" inparam_basic |awk '{print $2}' `

if ( $simtype == 'single' ) then
    set simdir1 = './'
    if ( `grep FINISHED OUTPUT_* | wc -l` == 0 ) then
        echo "SOLVER run not yet finished...hang on a second"
       exit
    endif
else if ( $simtype == 'moment' ) then
    set simdir1 = 'MZZ'
    if ( `grep FINISHED MZZ/OUTPUT_* | wc -l` == 0 ) then
        echo "SOLVER run MZZ not yet finished...hang on a second"
        exit
    endif
    if ( `grep FINISHED MXX_P_MYY/OUTPUT_* | wc -l` == 0 ) then
        echo "SOLVER run MXX_P_MYY not yet finished...hang on a second"
        exit
    endif
    if ( `grep FINISHED MXY_MXX_M_MYY/OUTPUT_* | wc -l` == 0 ) then
        echo "SOLVER run MXY_MXX_M_MYY not yet finished...hang on a second"
        exit
    endif
    if ( `grep FINISHED MXZ_MYZ/OUTPUT_* | wc -l` == 0 ) then
        echo "SOLVER run MXZ_MYZ not yet finished...hang on a second"
        exit
    endif
else
    echo 'postprocessing only implemented for SIMULATION_TYPE single and moment'
    exit
endif

set outdir = `grep "DATA_DIR" param_post_processing |awk '{print $2}' |sed 's/"/ /g' `
echo "All post-processed data in: "$outdir

if ( ! -d $outdir ) then
    mkdir $outdir
else
    echo 'ERROR: Output Directory already exists:'
    echo $outdir
    exit
endif

mkdir $outdir/SNAPS

cp -p param_* $outdir 
/bin/cp -p -f param_* $outdir
cp -p post_processing.csh $outdir
cp -p xpost_processing $outdir
cp -p post_processing.F90 $outdir 

cp -p plot_record_section.m $outdir

echo
echo "%%%%%%%%% PROCESSING seismograms/wavefields %%%%%%%%%"
mkdir $outdir/SEISMOGRAMS
mkdir $outdir/SEISMOGRAMS/UNPROCESSED

echo ".... output in "$outdir"/OUTPUT_postprocessing ...."

./xpost_processing > $outdir/OUTPUT_postprocessing

if ( $status != 0 ) then
    echo 'xpost_processing exited with error'
    echo 'check output for details:'
    echo $outdir'/OUTPUT_postprocessing'
    exit
endif

echo "Done with post processing, results in SEISMOGRAMS/ " 

if ( -f param_snaps) then 
    echo " .... and SNAPS/"
endif

set gnu_query = `which gnuplot | wc -l `
if ( $gnu_query == 1 ) then
    echo
    echo "%%%%%%%%% PLOTTING seismograms (gnuplot) %%%%%%%%%%"
    cd $outdir
    set seistype = "disp"
    echo "seismogram type:" $seistype
    set reclist = `cat $homedir/$simdir1/Data/receiver_names.dat |awk '{print $1}'`
    echo "1st receiver:" $reclist[1]
    set colat = `cat $homedir/$simdir1/Data/receiver_names.dat |awk '{print $2}' |sed 's/00000/ /g' |awk '{print $1}'`
    set lon = `cat $homedir/$simdir1/Data/receiver_names.dat |awk '{print $3}' |sed 's/00000/ /g'  |awk '{print $1}'`
    set epidist = `cat $homedir/$simdir1/Data/receiver_pts.dat |awk '{print $1}'`
    echo "1st receiver colatitude/longitude/epidist:" $colat[1] " " $lon[1] " " $epidist[1]

    set reccomp = `ls SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*.dat |sed 's/mij_/ /g ' |awk '{print $2}' | sed 's/_/ /g' |awk '{print $2}'  |sed 's/\.dat/ /g '`

    set conv = `ls SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*.dat |sed 's/_mij_/ /g '  | sed 's/\.dat/ /g ' |awk '{print $2}' `
    set t2 = `tail -n 1 SEISMOGRAMS/{$reclist[1]}_{$seistype}_post_mij_*{$reccomp[1]}.dat |awk '{print $1}' `
    echo "convolution:" $conv
    echo "receiver components:" $reccomp
    mkdir GRAPHICS
    set i = 0
    foreach  rec (${reclist})
    @ i++
    set j = 0
    foreach comp (${reccomp}) 
    @ j++
        set recname = `echo $rec"_"$seistype"_post_mij_"{$conv[$j]}`
        echo "Plotting receiver " $recname 
        echo 'set term png linewidth 1  ' >! plot_recs.plot
        echo 'set output "GRAPHICS/'$recname'.png"' >> plot_recs.plot
        echo 'set title "colat,lon: '$colat[$i], $lon[$i]', epidist: '$epidist[$i]'"'>> plot_recs.plot
        echo 'plot "SEISMOGRAMS/'$recname'.dat" with lines' >> plot_recs.plot
        echo "set xrange [ 0: "$t2"];set xlabel 'time [s]';set ylabel 'displacement [m]' " >> plot_recs.plot
        gnuplot plot_recs.plot
        cd GRAPHICS; 
        convert $recname.png $recname.gif
        convert $recname.png $recname.pdf
        rm -f $recname.png
        cd ..
    end
    end
    echo "Done with plotting, results in GRAPHICS/"
endif

cd $homedir
set taup_query = `which taup_time | wc -l `
if ( $taup_query == 1 ) then
    echo
    echo "%%%%%%%%% Computing TRAVELTIMES (taup) %%%%%%%%%"
    set depth = `grep "source depth" $simdir1/simulation.info |awk '{print $1}'`
    set model = `grep Background $simdir1/mesh_params.h |awk '{print $5}'`

    if ( $model == 'prem' || $model == 'iasp91' ) then
	set num_rec = `wc -l $simdir1/Data/receiver_pts.dat |awk '{print $1}'`
	set epi_list = `tail -n $num_rec $simdir1/Data/receiver_pts.dat | awk '{print $1}'`
	set depth_short = `echo $depth |sed 's/\./ /g' |awk '{print $1}'` 
	echo "Earthquake depth:" $depth_short
	set i = 0
	cd $outdir
	foreach rec (${epi_list})
	    @ i++
	    echo "traveltimes for epicentral distance" $rec

	    set tt =  `taup_time -mod $model -h $depth -ph pP -deg $rec | grep " pP " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_pP_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph PP -deg $rec | grep " PP " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_PP_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph SS -deg $rec | grep " SS " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_SS_traveltime2.dat ;	    endif

	    if ( $rec < 100. ) then
	    set tt =  `taup_time -mod $model -h $depth -ph P -deg $rec | grep " P " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_P_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph S -deg $rec | grep " S " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_S_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph PcP -deg $rec | grep " PcP " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_PcP_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph ScS -deg $rec | grep " ScS " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_ScS_traveltime2.dat ;	    endif
	    endif

	    if ( $rec > 95. ) then
	    set tt =  `taup_time -mod $model -h $depth -ph Pdiff -deg $rec | grep " Pdiff " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_Pdiff_traveltime2.dat ;	    endif

	    set tt =  `taup_time -mod $model -h $depth -ph Sdiff -deg $rec | grep " Sdiff " | awk '{print $4}' |head -n 1 |grep -v "=="` 
	    if ( ${#tt} == 1) then ;	    echo $tt $rec >> taup_Sdiff_traveltime2.dat ;	    endif

	    endif

	end
	sort -n taup_P_traveltime2.dat |grep -v "==>" > taup_P_traveltime.dat
	sort -n taup_pP_traveltime2.dat |grep -v "==>" > taup_pP_traveltime.dat
	sort -n taup_S_traveltime2.dat |grep -v "==>" > taup_S_traveltime.dat
	sort -n taup_PP_traveltime2.dat |grep -v "==>" > taup_PP_traveltime.dat
	sort -n taup_SS_traveltime2.dat |grep -v "==>" > taup_SS_traveltime.dat
	sort -n taup_PcP_traveltime2.dat |grep -v "==>" > taup_PcP_traveltime.dat
	sort -n taup_ScS_traveltime2.dat |grep -v "==>" > taup_ScS_traveltime.dat
	sort -n taup_Pdiff_traveltime2.dat |grep -v "==>" > taup_Pdiff_traveltime.dat
	sort -n taup_Sdiff_traveltime2.dat |grep -v "==>" > taup_Sdiff_traveltime.dat
	mkdir TAUP
	rm -f taup_*2.dat
	mv taup_*.dat TAUP
	echo "Done with taup, results in TAUP/"
    endif
endif

echo
echo
echo " %%%%%%%%%%% QUICK OVERVIEW OF RESULTS %%%%%%%%%%%"
echo "1) Check earthquake and seismograms: "
echo "        open google earth"
echo "        load file "$outdir"/googleearth_src_rec_seis.kml"
echo "        locate earthquake, click and check parameters"
echo "        click on receiver location, check location and seismograms"
echo
echo " 2) Seismogram record section:"
echo "        open matlab in directory" $outdir
echo "        >> plot_record_section"
echo "        plots of seismograms sorted in epicentral distance with relative amplitudes"
echo "        station name given on the left of each trace"
echo "        traveltimes for some phases (if simulations are upon PREM or IASP91)"
echo
echo " 3) Seismogram time series:"
echo "        load individual seismograms, e.g. with the command:"
echo "        >> xmgrace" $outdir"/SEISMOGRAMS/"$recname".dat"

if ( -f param_snaps ) then
    echo
    echo " 4) Wavefield snapshots:"
    echo "        open paraview"
    echo "        load "$outdir"/SNAPS/snap_mij_cell_*vtk"
    echo "        ... this should be a wavefield movie in 3D. "
    echo "        adjust normalization"
endif

echo
echo "==========================================="
echo "                  DONE with post processing."
echo "==========================================="
