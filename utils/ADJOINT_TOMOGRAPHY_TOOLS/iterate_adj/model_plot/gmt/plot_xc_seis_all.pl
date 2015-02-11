#!/usr/bin/perl -w

#==========================================================
#
#  plot_xc_seis_all.pl
#  Carl Tape
#  27-March-2009
#
#  EXAMPLES:
#     plot_xc_seis_all.pl T006_T030 1 1
#     plot_xc_seis_all.pl T006_T030 0 1
#     plot_xc_seis_all.pl T003_T030 1 1
#     plot_xc_seis_all.pl T003_T030 0 1
#     plot_xc_seis_all.pl T002_T030 1 1
#     plot_xc_seis_all.pl T002_T030 0 1
#
#==========================================================

if (@ARGV < 3) {die("Usage: plot_xc_seis_all.pl xxx\n");}
($Ttag,$iwin,$i1D) = @ARGV;

# make output directory
$otag = "${Ttag}_win${iwin}_1D${i1D}";
$odir = "SAVE/$otag";
#`rm -rf $odir`;
`mkdir -p $odir`;

# directory containing the data files
$irun = 1;
$stirun = sprintf("%2.2i",$irun);

$dir0 = "/home/carltape/ADJOINT_TOMO";
$dirin = "${dir0}/ADJOINT_TOMO_OUTPUT";
$dirdat = "${dirin}/model_plot_matlab_OUTPUT/vert_${stirun}";
$seisdir = "${dirin}/misfit_plot_OUTPUT/${otag}";
$xcdir = "${dirin}/model_plot_gmt_OUTPUT/VERT_${stirun}";
if (not -e $dirdat) {die("Check if dirdat $dirdat exist or not\n");}
if (not -e $xcdir) {die("Check if xcdir $xcdir exist or not\n");}
if (not -e $seisdir) {die("Check if seisdir $seisdir exist or not\n");}

#---------------------------------------------------------

# load fault positions along each profile
$ffile = "$dirdat/ALL_rays_fault_positions_mod";
if (not -f $ffile) {die("Check if ffile $ffile exist or not\n");}
open(IN,$ffile); @faults = <IN>; close(IN);
$nump = @faults;

for ($p = 1; $p <= $nump; $p ++ ) {
   (undef,$eid,$stanet,$azi1,$azi2,$saf,$gf,$sgf,$mcf,$syf,$crf,$ef,$snf,$kcf) = split(" ",$faults[$p-1]);
   print "-- $p -- $stanet -- $eid --\n";
}
#die("TESTING");

$pmin = 1; $pmax = $nump;
#$pmin = 1; $pmax = 10;
$pmin = 77; $pmax = $pmin;

for ($p = $pmin; $p <= $pmax; $p ++ ) {
   $stp = sprintf("%3.3i",$p);

   # positions of six faults along each profile
   # (1) SAF, (2) GF, (3) SGF, (4) MCF, (5) SYF, (6) CRF, (7) EF, (8) SN, (9) KC
   (undef,$eid,$stanet,$azi1,$azi2,$saf,$gf,$sgf,$mcf,$syf,$crf,$ef,$snf,$kcf) = split(" ",$faults[$p-1]);
   ($sta,$net) = split("\\.",$stanet);    # only for irun = 1
   $sazi1 = sprintf("%3.3i",$azi1);
   $sazi2 = sprintf("%3.3i",$azi2);
   print "-- $p -- $stanet -- $sta -- $net -- $eid -- $sazi1 -- $sazi2 --\n";

   # length of cross section
   $rayfile = "$dirdat/vert_xc_${stp}_ray_path";
   if (not -f $rayfile) {die("Check if rayfile $rayfile exist or not\n");}
  (undef,undef,undef,undef,undef,undef,undef,undef,$dmin,$dmax) = split(" ",`minmax -C -I2 $rayfile`);
  $dran = ($dmax - $dmin)/1000;
  $heightxc = $heightxc0;
  if ($dran < 350) {$heightxc = "width=4.5cm";}
  else {$heightxc = "height=12cm";}
  print "length of cross section is $dran km\n";

   # find cross section file
   @xcfiles = glob("${xcdir}/*$sta*$net*$eid*eps");
   $nxcfiles = @xcfiles;
   if($nxcfiles != 1) {
      print "$ftag\n @xcfiles\n $nxcfiles\n";
      die("Must have only one cross section file");
   }
   $xcfile = $xcfiles[0];
   #$xcfile = "${xcdir}/vert_01_xc_LDF_CI_269_9983429_vs_m16_m00_001.eps";
   if (not -f $xcfile) {die("Check if xcfile $xcfile exist or not\n");}

   # find seismofile
   $stag = "seis_2";
   if($i1D==1) {$stag = "seis_3";}
   $ftag = "${seisdir}/${eid}_${Ttag}_${sta}_${net}*${stag}.eps";
   @seisfiles = glob($ftag);
   $nseisfiles = @seisfiles;
   if($nseisfiles != 1) {
      print "$ftag\n @seisfiles\n $nseisfiles\n";
      die("Must have only one seismogram file");
   }
   $seisfile = $seisfiles[0];
   #$seisfile = "${seisdir}/9983429_T006_T030_LDF_CI_m16_seis_2.eps";
   if (not -f $seisfile) {die("Check if seisfile $seisfile exist or not\n");}

   # write latex file for combining figures
   $name = "xc_and_seis";
   $texfile = "${name}.tex";
   open(TEX,">$texfile");
print TEX "\\documentclass[pdf,mpa]{prosper}\n";
print TEX "\\usepackage{portland}\n";
print TEX "\\begin{document}\n";
print TEX "\\begin{slide}{}\n";
print TEX "\\vspace{-1.7cm}\n";
print TEX "\\fcolorbox{white}{white}{\n";
print TEX "\\hspace{-1cm} \n";
print TEX "\\begin{tabular}{c}\n";
#print TEX "\\includegraphics[width=3.8cm,angle=-90]{$xcfile}\n";
print TEX "\\includegraphics[${heightxc},angle=-90]{$xcfile}\n";
print TEX "\\\\\n";
print TEX "\\includegraphics[width=3.8cm,angle=-90]{$seisfile}\n";
print TEX "\\end{tabular}\n";
print TEX "}\n";
print TEX "\\end{slide}\n";
print TEX "\\end{document}\n";
close(TEX);

   # execute latex file and convert to pdf
   $fname = "xc_and_seis_${eid}_${sta}_${net}_${Ttag}";
   `latex ${name}`;
   `dvips -o ${name}.ps ${name}`;
   `mv ${name}.ps ${fname}.ps`;
   `ps2pdf ${fname}.ps`;

}

# UNCOMMENT: move all pdf and ps files to output directory
#`mv xc_and_seis*pdf xc_and_seis*ps $odir`;
#`/home/carltape/bin/pdcat -r $odir/xc_and_seis*pdf $odir/ALL_xc_and_seis_${otag}.pdf`;

#==================================================
