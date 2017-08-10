#!/bin/csh -f

set list_grid_files = `ls messagepassing_solid_1_*.dat`
foreach gridfile (${list_grid_files})
set procnum = `ls $gridfile | sed 's/_1_/ /' | sed 's/\.dat/ /' | awk '{print $2}' `
psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'MPI solid f=1':WeSn -U/0/-1/ALX"`pwd`" -K -V >mpitest_sol1_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> mpitest_sol1_proc$procnum.ps
echo wrote processor valence into mpitest_sol1_proc$procnum.ps
end

set list_grid_files = `ls messagepassing_solid_2_*.dat`
foreach gridfile (${list_grid_files})
set procnum = `ls $gridfile | sed 's/_2_/ /' | sed 's/\.dat/ /' | awk '{print $2}' `
psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'MPI solid f=2':WeSn -U/0/-1/ALX"`pwd`" -K -V >mpitest_sol2_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> mpitest_sol2_proc$procnum.ps
echo wrote processor valence into mpitest_sol2_proc$procnum.ps
end

set list_grid_files = `ls messagepassing_solid_3_*.dat`
foreach gridfile (${list_grid_files})
set procnum = `ls $gridfile | sed 's/_3_/ /' | sed 's/\.dat/ /' | awk '{print $2}' `
psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'MPI solid f=3':WeSn -U/0/-1/ALX"`pwd`" -K -V >mpitest_sol3_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> mpitest_sol3_proc$procnum.ps
echo wrote processor valence into mpitest_sol3_proc$procnum.ps
end

set list_grid_files = `ls messagepassing_fluid_*.dat`
foreach gridfile (${list_grid_files})
set procnum = `ls $gridfile | sed 's/_fluid_/ /' | sed 's/\.dat/ /' | awk '{print $2}' `
psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'MPI fluid f=3':WeSn -U/0/-1/ALX"`pwd`" -K -V >mpitest_fluid_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> mpitest_fluid_proc$procnum.ps
echo wrote processor valence into mpitest_fluid_proc$procnum.ps
end

