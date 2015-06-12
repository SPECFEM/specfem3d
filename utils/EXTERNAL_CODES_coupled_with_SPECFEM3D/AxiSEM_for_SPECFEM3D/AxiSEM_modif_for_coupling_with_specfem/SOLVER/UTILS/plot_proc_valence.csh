#!/bin/csh -f

####### VALENCE FOR SOLID DOMAIN PER PROCESSOR #############
set list_grid_files = `ls valence_solid_*.dat`
foreach gridfile (${list_grid_files})
set procnum = `ls $gridfile | sed 's/solid_/ /' | sed 's/\.dat/ /' | awk '{print $2}' `
psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'Valence solid':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_solid_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> valence_solid_proc$procnum.ps
echo wrote processor valence into valence_solid_proc$procnum.ps
end

####### VALENCE FOR FLUID DOMAIN PER PROCESSOR #############
set list_grid_files = `ls valence_fluid_*.dat`
foreach gridfile (${list_grid_files})
set procnum = `ls $gridfile | sed 's/fluid_/ /' | sed 's/\.dat/ /' | awk '{print $2}' `
psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'Valence fluid':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_fluid_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> valence_fluid_proc$procnum.ps
echo wrote processor valence into valence_fluid_proc$procnum.ps
end
