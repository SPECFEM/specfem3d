psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence.dat -P \
-B::/:::.'Global valence ':WeSn -U/0/-1/"`pwd`" -K -V >valence_glob.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_glob.ps

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence_fluid.dat -P \
-B::/:::.'Fluid valence ':WeSn -U/0/-1/"`pwd`" -K -V >valence_fluid.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_fluid.ps

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence_solid.dat -P \
-B::/:::.'Solid valence':WeSn -U/0/-1/"`pwd`" -K -V >valence_solid.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_solid.ps

