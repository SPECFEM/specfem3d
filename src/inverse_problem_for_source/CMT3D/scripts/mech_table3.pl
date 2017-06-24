#!/usr/bin/perl
use lib '/opt/seismo/lib/perl';
use GMT_PLOT;
use GMT_PLOT_SC;
use CMT_TOOLS;
use POSIX;
use List::Util qw[max];

## ok, this version includes both 2, 3 and 6 seconds measurements
# parameter files
$inv_file="INVERSION.PAR.SAVE";
$outdir="mech"; $psfile="mech.ps";

$weigh_data="2 2 1 0.5 1.15 0.55 0.78";
$weigh_data_nocomp="1 1 1 0.5 1.15 0.55 0.78";
$weigh_data_noaz="2 2 1 0 1.15 0.55 0.78";

if (not -f "CMTSOLUTION") {die("Check cmt solution file\n");}
($evid)=split(" ",`cat CMTSOLUTION | awk 'NR == 2 {print \$3}'`);
$mwin="MEASUREMENT_WINDOWS";

@names=("Initial","Grid Search","7 Par+ZT",
        "7 Par+no-SC","6 par+DC","7 Par+DC",
        "7 Par(2-30s)","7 Par(3-30s)","7 Par(6-30s)");
@npars=(6,6,7, 7,6,7, 7,7,7);
@wdatas=($weigh_data, $weigh_data, $weigh_data, $weigh_data,$weigh_data,$weigh_data, $weigh_data,$weigh_data,$weigh_data);
@scorrs=("true","true",".true.",  ".false.", ".true.", ".true.",  ".true.",".true.",".true.");
@cons=(".true. true. 0.0", ".true. .true. 0.0", ".true. .false. 0.0",
       ".true. .false. 0.0", "true. .false. 0.0", "true. .true. 0.0",
       ".true. .false. 0.0",".true. .false. 0.0",".true. .false. 0.0");
@meas=("${mwin}_${evid}_ALL","${mwin}_${evid}_ALL","${mwin}_${evid}_ALL",
       "${mwin}_${evid}_ALL","${mwin}_${evid}_ALL","${mwin}_${evid}_ALL",
       "${mwin}_${evid}_T002_T030_m12","${mwin}_${evid}_T003_T030_m12","${mwin}_${evid}_T006_T030_m12");

$nbeach=@npars;
$ncols = 3;
$nrows = 3;
# plot size jx/jy, plot origin in xy
($jx,$jy) = shift_xy("1 1","$ncols $nrows","7  7", \@xy,"0.8 0.8");
$JX = "-JX$jx/$jy"; $R="-R0/3/0/3";
$GMT_PLOT::paper_orient="-P";

if (not -f "cmt3d_flexwin") {die("Check if cmt3d_flexwin exists or not\n");}

open(INV,"$inv_file")|| die("check $inv_file\n");
@all=<INV>;close(INV);

$cmt=$all[0]; $cmt_new=$all[1]; $npar=$all[2]; $delta=$all[3];
$flex=$all[4];$lwdata=$all[5]; $wdata=$all[6]; $scorr=$all[7]; $con=$all[8];
$wns=$all[9];
chomp($cmt_new);chomp($cmt);
(undef,undef,$ename)=split(" ",`grep 'event name' $cmt`);

@cmt=($cmt);@cmt=(@cmt,"${cmt}_GRD");
if (-f "${cmt}_GRD") {print BASH "cp -f ${cmt}_GRD $outdir\n";}
for ($i=2;$i<$nbeach;$i++) {@cmt=(@cmt,"CMTSOLUTION.$i");}

open(BASH,">mech_table.bash");
print BASH "mkdir -p $outdir\n";
print BASH "cp -f $cmt $outdir\n";
if (-f "${cmt}_GRD") {print BASH "cp -f ${cmt}_GRD $outdir\n";}

for ($k=2;$k<$nbeach;$k++) { # cmt3d runs
  print BASH "echo Running cmt3d_flexwin for INVERSION.$k ...\n";
  open(INV,">INVERSION.$k");
  print INV "${cmt}\n${cmt_new}\n$npars[$k]\n$delta$meas[$k]\n.true.\n$wdatas[$k]\n$scorrs[$k]\n$cons[$k]\n";
  if ($k==2) {print INV ".true.\n";}   else { print INV "${wns}\n";}
  close(INV);
  print BASH "# --- mech $k ----\n";
  print BASH "cp -f INVERSION.$k INVERSION.PAR\n";
  print BASH "cmt3d_flexwin > cmt3d.stdout\n";
  print BASH "if [ \$? != 0 ]; then\n  echo \"exiting case $k\";exit\nfi\n";
  print BASH "perl -pi -e \"s/INV/$names[$k]/\" $cmt_new\n";
  print BASH "mv -f $cmt_new $outdir/$cmt[$k]\n";
  print BASH "mv -f cmt3d_flexwin.out $outdir/cmt3d_out.$k\n";
  print BASH "mv -f cmt3d.stdout $outdir/cmt3d_stdout.$k\n";
  print BASH "mv -f INVERSION.$k $outdir\n";
}
close(BASH);
system("chmod a+x mech_table.bash; mech_table.bash");
for ($k=2;$k<$nbeach;$k++) {
  (@tmp)=split(" ",`grep Variance $outdir/cmt3d_stdout.$k`);
  $var[$k]=$tmp[8]; print "$k -- VR = $var[$k]\n";}

#die("here\n");
# get depth information to display
($odep)=split(" ",`cmtsol2faultpar.pl $cmt | grep mw | cut -d = -f 2 | cut -d / -f 2`);
($ndep)=split(" ",`cmtsol2faultpar.pl $outdir/$cmt[2] | grep mw | cut -d = -f 2 | cut -d / -f 2`);
($dep1)=split(" ",`cmtsol2faultpar.pl $outdir/$cmt[6] | grep mw | cut -d = -f 2 | cut -d / -f 2`);
($dep2)=split(" ",`cmtsol2faultpar.pl $outdir/$cmt[7] | grep mw | cut -d = -f 2 | cut -d / -f 2`);
($dep3)=split(" ",`cmtsol2faultpar.pl $outdir/$cmt[8] | grep mw | cut -d = -f 2 | cut -d / -f 2`);
$ddep=sprintf("%5.2f",max(abs($dep1-$ndep),abs($dep2-$ndep),abs($dep3-$ndep)));

open(BASH,">mech_plot.bash");
print BASH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 10 MEASURE_UNIT inch PAPER_MEDIA letter TICK_LENGTH 0.1c\n";

plot_psxy(\*BASH,$psfile,"$JX -K -X0 -Y0 $R","");
$k=0;
for ($i=0;$i<$nrows;$i++) {
  for ($j=0;$j<$ncols;$j++) {
    print BASH "# ---- Mech $k, $cmt[$k]------\n";
    print "---- Mech $k, $cmt[$k]------\n";
    $xy=$xy[$i][$j];($x,$y) = split(/\//,$xy);
    plot_psxy(\*BASH,$psfile,"$JX -X$x -Y$y $R","");
    if (-f "$outdir/$cmt[$k]") {
      plot_psmeca_raw(\*BASH,$psfile,"-JX -R -Sm1.0 ","1.5 2","$outdir/$cmt[$k]");
      @output = `cmtsol2faultpar.pl $outdir/$cmt[$k]`;
      (undef,$tmp1) = split(" ", $output[5]);
      (undef,$tmp2) = split(" ", $output[7]);
      ($s,$d,$r,$p,$t) = split(/\//,$tmp2);
      $tex="1 1 9 0 4 CM $names[$k]\n";
      $tex.="0.2 0.8 9 0 4 LM Mw/dep/eps/= $tmp1\n0.2 0.6 9 0 4 LM S/D/R= $s/$d/$r";
      if ($k>=2) {$tex.=sprintf("\n0.2 0.4 9 0 4 LM Var=%6.2f%",$var[$k]);}
      plot_pstext(\*BASH,$psfile,"-JX -R  -G0/0/255","$tex");
    }
    if ($k==1) {
     ($nt) = split(" ",`grep 'T.sac.d' $meas[2] | wc | awk '{print \$1}'`);
     ($nr) = split(" ",`grep 'R.sac.d' $meas[2] | wc | awk '{print \$1}'`);
     ($nz) = split(" ",`grep 'Z.sac.d' $meas[2] | wc | awk '{print \$1}'`);
     ($n2) = split(" ",`grep 'T002.*.sac.d' $meas[2] | wc | awk '{print \$1}'`);
     ($n3)  = split(" ",`grep 'T003.*.sac.d' $meas[2] | wc | awk '{print \$1}'`);
     ($n6)  = split(" ",`grep 'T006.*.sac.d' $meas[2] | wc | awk '{print \$1}'`);

      plot_pstext(\*BASH,$psfile,"-JX -R  -N","0.5 4 12 0 4 CM $ename ($nz/$nr/$nt; $n2/$n3/$n6; $odep/$ndep/$ddep)");}
    plot_psxy(\*BASH,$psfile,"-JX -X-$x -Y-$y","");
    $k++;
  }
}
plot_psxy(\*BASH,$psfile,"-JX -O","");
close(BASH);
system("chmod a+x mech_plot.bash; mech_plot.bash");
