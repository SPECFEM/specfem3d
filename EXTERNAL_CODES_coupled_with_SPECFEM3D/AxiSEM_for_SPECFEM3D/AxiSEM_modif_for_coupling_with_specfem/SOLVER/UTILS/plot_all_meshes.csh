#!/bin/csh -f

########### SOLID MESHES PER PROCESSOR ####################
set list_grid_files_solid = `ls serend_mesh_solid.dat0*`

foreach gridfile (${list_grid_files_solid})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_solid_proc{$procnum}.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"$gridfile\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
echo wrote processor grid to grid_solid_proc{$procnum}.ps

end


########### FLUID MESHES PER PROCESSOR ####################
set list_grid_files_fluid = `ls serend_mesh_fluid.dat0*`

foreach gridfile (${list_grid_files_fluid})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_fluid_proc{$procnum}.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"$gridfile\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
echo wrote processor grid to grid_fluid_proc{$procnum}.ps

end

########### GLOBAL SOLID MESH ####################

if ( -f serend_mesh_solid.dat ) then
rm -f serend_mesh_solid.dat
endif

cat serend_mesh_solid*dat0*  > serend_mesh_solid.dat

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_solid_solver.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"serend_mesh_solid.dat\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
rm -f serend_mesh_solid.dat
echo wrote processor grid to grid_solid_solver.ps


########### GLOBAL FLUID MESH ####################

if ( -f serend_mesh_fluid.dat ) then
rm -f serend_mesh_fluid.dat
endif

cat serend_mesh_fluid*dat0*  > serend_mesh_fluid.dat

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_fluid_solver.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"serend_mesh_fluid.dat\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
rm -f serend_mesh_fluid.dat
echo wrote processor grid to grid_fluid_solver.ps




########### GLOBAL JOINT MESH ####################

if ( -f serend_mesh.dat ) then
rm -f serend_mesh.dat
endif

cat serend_mesh_*dat*  > serend_mesh.dat

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_solver.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"serend_mesh.dat\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
rm -f serend_mesh.dat
echo wrote processor grid to grid_solver.ps






