psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence_meshtest.dat -P \
-B::/:::.'Valence meshtest':WeSn -U/0/-1/"`pwd`" -K -V >valence_meshtest.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_meshtest.ps
#gs valence.ps
