#!/bin/csh -f

set list_grid_files_solid = `ls serend_mesh_solid.dat0*`

foreach gridfile (${list_grid_files_solid})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_proc{$procnum}.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"$gridfile\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
echo wrote processor grid to grid_solid_proc{$procnum}.ps

end

set list_grid_files_fluid = `ls serend_mesh_fluid.dat0*`

foreach gridfile (${list_grid_files_fluid})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

if ( -f plot_oneproc_mesh.plot) then
rm -f plot_oneproc_mesh.plot
endif

echo 'set size ratio -1' > plot_oneproc_mesh.plot
echo 'set term postscript portrait solid "Helvetica" 24' >> plot_oneproc_mesh.plot
echo set output \"grid_proc{$procnum}.ps\" >> plot_oneproc_mesh.plot
echo 'set noborder'>> plot_oneproc_mesh.plot
echo 'set noxtics ; set noytics'>> plot_oneproc_mesh.plot
echo plot \"$gridfile\" t\'\' with l lw .4 >> plot_oneproc_mesh.plot

gnuplot plot_oneproc_mesh.plot
echo wrote processor grid to grid_fluid_proc{$procnum}.ps

end





