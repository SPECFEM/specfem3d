~/SVN/mt_measure_adj/scripts/mtm.pl -s SYNTH,qmxd.sac.bp50.500 -o OUTPUT_FILES/ -w t3/300/300 -V DATA/2001.026.0*bp50.500

~/SVN/mt_measure_adj/scripts/mtm.pl -s SYNTH,qmxd.sac.bp50.500 -o OUTPUT_FILES/ -w t3/300/300 -V -C DATA/2001.026.0*bp50.500

~/SVN/mt_measure_adj/scripts/plot_mtm.pl OUTPUT_FILES/*.mtm.dlnA

~/SVN/mt_measure_adj/scripts/plot_adj.pl OUTPUT_FILES/*.mtm.adj
