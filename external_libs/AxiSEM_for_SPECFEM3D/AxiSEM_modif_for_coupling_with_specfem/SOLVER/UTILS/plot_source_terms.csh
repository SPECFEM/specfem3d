psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 src_term_norm1.dat -P \
-B::/:::.'source term 1':WeSn -U/0/-1/"`pwd`" -K -V >source_term1.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>source_term1.ps

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 src_term_norm2.dat -P \
-B::/:::.'source term 2':WeSn -U/0/-1/"`pwd`" -K -V >source_term2.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>source_term2.ps

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 src_term_norm3.dat -P \
-B::/:::.'source term 3':WeSn -U/0/-1/"`pwd`" -K -V >source_term3.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>source_term3.ps

