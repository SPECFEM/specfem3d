#!/usr/bin/perl -w

# this program plots sophisticated data and syn waveforms according to either
# given data and syn dir
# or flexwin output file 
# sort by distance and then by azimuth, a regional SC map is also supplied

use lib '/home/liuqy/local/perl-lib';
use GMT_PLOT;
use GMT_PLOT_SC;
#use CMT_TOOLS;
use POSIX;
use Getopt::Std;
use List::Util qw[min max];


if (@ARGV == 0) {die("Usage: plot_data_and_syn.pl [-d data_dir,data_ext -s syn_dir,syn_ext -Snew_syn_dir,new_syn_ext]/-M flexwin_file -m CMTSOLUTION -A Z/R/T(0.03/0.03/0.01) -C\n");}

if (!getopts('d:s:S:m:M:A:C')) {die("Check input options\n");}

# data/syn dir and extension
if ($opt_M) {
  $flexwin=1; $mfile = $opt_M;
} elsif ($opt_d or $opt_s) {
  $flexwin=0;
  if ($opt_d) {
    ($data_dir,$data_ext) = split(/\,/,$opt_d);
    if ($data_ext) {$data_ext = ".$data_ext";}
    else {$data_ext="";}}
  if ($opt_s) {
    ($syn_dir,$syn_ext) = split(/\,/,$opt_s);
    if ($syn_ext) {$syn_ext=".$syn_ext";}
    else {$syn_ext="";}}
} else {die("Input -M or -d/-s options\n");}

if ($opt_S) {$plot_new = 1;
  ($new_syn_dir,$new_syn_ext) = split(/\,/,$opt_S);
  if ($new_syn_ext) {$new_syn_ext=".$new_syn_ext";}
  else {$new_syn_ext="";}}

if ($opt_m) {$cmt_file=$opt_m;} else {die("Input cmt file \n");}
($moment) = split(" ",`cmtsol2faultpar.pl $cmt_file | awk 'NR == 3 {print \$3}'`);

if ($opt_c) {$socal = 1;} else {$socal = 0;}
# comps
$ncols=3; @comps=("Z","R","T");#Z R T

# this controls the paper orientation;
$GMT_PLOT::paper_orient = "-P";

# *********** Modify: CONTROL PARAMETERS *************
$nfiles=12; # this is number of files per plot

($jx,$jy) = shift_xy("1 1","$ncols 1","7  7", \@xy,"0.95 0.95");
print "xy = $xy[0][0], $xy[0][1], $xy[0][2]\n";
$JX = "-JX$jx/$jy"; 
if ($socal) {$B1="-B20/4";} else {$B1="-B1000/4";}
@B=("${B1}WSN","${B1}SN","${B1}ESN");

$Vl=3.0; # surface wave velocity to determine tend
$bp=10; $bp2=10;# start time before t1 header and end time after tend
$E="-Ent1"; # align by trace number and t1
if ($opt_A) {
  ($za,$ra,$ta) = split(/\//,$opt_A);
  @size_all=($za/2.4e23*$moment,$ra/2.4e23*$moment,$ta/2.4e23*$moment);
} else {
  @size_all=(0.03/2.4e23*$moment,0.03/2.4e23*$moment,0.06/2.4e23*$moment);}
print "plot sizes = @size_all \n";
$wdata="-W3/0/0/0"; $wsyn="-W2/255/0/0";  $gwin="-G220 -W3";
if ($plot_new) {$wsyn="-W2/0/0/255"; $wnsyn="-W2/255/0/0";}

$jxy_map="-X3 -Y8"; $JM="-JM3";
# *****************************************************

# loop over data files, select corresponding syn files and set up % file
print "Matching data files with synthetics\n";%file=();
if ($flexwin == 0) { # given data and syn dir
  @files = glob("${data_dir}/*H?$data_ext");

  foreach $datfile (@files) {
    (undef,$net,$sta,$chan,$dist,$az,$t1,$stla,$stlo) = split(" ",`saclst knetwk kstnm kcmpnm dist az t1 stla stlo f $datfile`);
    ($comp) = split(" ",`echo $chan | awk '{print substr(\$1,3,1)}'`);
    $synfile = "${syn_dir}/${sta}.${net}.BH${comp}${syn_ext}";
    
    if (not -f $synfile) {
      print "**** Missing synthetic pair $synfile for $datfile ***, skipping ...\n"; } 
    else {
      $file{$sta}{dist}=$dist; $file{$sta}{az}=$az; $file{$sta}{P}=$t1;
      $file{$sta}{"data-$comp"}=$datfile; $file{$sta}{"syn-$comp"}=$synfile;
      $file{$sta}{stla}=$stla; $file{$sta}{stlo}=$stlo;
      if ($plot_new) {
	$new_synfile="${new_syn_dir}/${sta}.${net}.BH${comp}${new_syn_ext}";
	$file{$sta}{"new-syn-$comp"}=$new_synfile;}
    }
  }
} else { # given flexwin output
  open(MFILE,"$mfile");
  $nfile=<MFILE>;
  for ($i=0;$i<$nfile;$i++) {
    # data and syn files
    ($datfile)=split(" ",<MFILE>); ($synfile)=split(" ",<MFILE>);
    if (not -f $datfile or not -f $synfile){die("Check $datfile or $synfile\n");}
    (undef,$net,$sta,$chan,$dist,$az,$t1,$stla,$stlo) = split(" ",`saclst knetwk kstnm kcmpnm dist az t1 stla stlo f $datfile`);
    ($comp) = split(" ",`echo $chan | awk '{print substr(\$1,3,1)}'`);
    (undef,$sta2,$chan2) = split(" ",`saclst kstnm kcmpnm f $datfile`);
    ($comp2) = split(" ",`echo $chan2 | awk '{print substr(\$1,3,1)}'`);
    if ($sta ne $sta2 or $comp ne $comp2) {die("Check data and syn pair $datfile and $synfile\n");}
    $file{$sta}{dist}=$dist; $file{$sta}{az}=$az; $file{$sta}{P}=$t1;
    $file{$sta}{"data-$comp"}=$datfile; $file{$sta}{"syn-$comp"}=$synfile;
    $file{$sta}{stla}=$stla; $file{$sta}{stlo}=$stlo;
    if ($plot_new) {
      $new_synfile="${new_syn_dir}/${sta}.${net}.BH${comp}${new_syn_ext}";
      $file{$sta}{"new-syn-$comp"}=$new_synfile;}
    # windows
    ($nwin)=split(" ",<MFILE>);
    $file{$sta}{"nwin-$comp"}=$nwin; $file{$sta}{"win-$comp"}="";
    for ($j=0;$j<$nwin;$j++) {
      $tse=<MFILE>; $file{$sta}{"win-$comp"}.=$tse;}
#    print "$datfile,$synfile,$sta,$nwin===\n";
  }
  close(MFILE);
}

# set up files in each plot
@sorted_stas=sort {$file{$a}{dist} <=> $file{$b}{dist}} keys %file;
$j=0;
for ($i=0;$i<@sorted_stas;$i++) {
  if ($i == ($j+1) * $nfiles) {$j++;}
  $plot[$j] .="$sorted_stas[$i] ";
}
$nplot=$j+1;

# loop over plots, plot data/syn files, window infomration, and maps
open(BASH,">plot_ds.bash");
print BASH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 10 MEASURE_UNIT inch PAPER_MEDIA letter TICK_LENGTH 0.1c\n";

for ($j=0;$j<$nplot;$j++) {
  # sort by dist
  @sta_dist=split(" ",$plot[$j]);
  # sort by azimuth
  @stas=sort {$file{$a}{az} <=> $file{$b}{az}} @sta_dist;
  @stlas=sort {$file{$a}{stla} <=> $file{$b}{stla}} @sta_dist;
  @stlos=sort {$file{$a}{stlo} <=> $file{$b}{stlo}} @sta_dist;

  # min dist and max dist
  $dist_min=$file{$sta_dist[0]}{dist}; $dist_max=$file{$sta_dist[-1]}{dist};

  # **** CONTROL PAR ***** (tmin/tmax, rel. to t1)
  $tmin=-1.5*$bp;
  $tmax=${dist_max}/$Vl-$file{$sta_dist[0]}{P}+2.5*$bp2;
  if ($tmax < 0) {die("Check $plot[$j] for tmin and tmax\n");}
  $C="-C$tmin/$tmax"; # cutting
  $tt1=$tmin-$bp/2; $tt2=$tmax+$bp2/2;
  $R="-R$tt1/$tt2/-1/$nfiles"; # region
  print "*** plot $j === @stas === $tt1/$tt2 === \n";

  print BASH "# ---- Plot $j ------\n";
  $psfile=sprintf("ds_%02d.ps",$j);
  plot_psxy(\*BASH,$psfile,"$JX $R -K -X0 -Y0","");

  for ($i=0;$i<$ncols;$i++) {
    $xy = $xy[0][$i];($x,$y) = split(/\//,$xy);
    plot_psxy(\*BASH,$psfile,"$JX -X$x -Y$y","");

    # ***** CONTROL PAR ******
    $size=2*$size_all[$i]/($dist_min+$dist_max);

    $data_files="";$syn_files="";$win_files=""; $new_syn_files="";
    if ($i==0) {$tex_file="";$xy_file="";$t_file="";}
    for ($k=0;$k<@stas;$k++) {
      $dfile=$file{$stas[$k]}{"data-$comps[$i]"};
      $sfile=$file{$stas[$k]}{"syn-$comps[$i]"};
      if ($plot_new) {$nsfile=$file{$stas[$k]}{"new-syn-$comps[$i]"};}
      # station name, dist, az
      if ($i==0) {
	$kt1=$k+0.2; $tt3=$tmin; $kt2 = $k-0.2;
	$dist=sprintf("%5.2f",$file{$stas[$k]}{dist});
	$az=sprintf("%5.2f",$file{$stas[$k]}{az});
	$tex_file.="$tt3 $kt1 8 0 4 LM $stas[$k] $dist\n$tt3 $kt2 8 0 4 LM $az\n";
	$stla=$file{$stas[$k]}{stla}; $stlo=$file{$stas[$k]}{stlo};
#	$stla2=$stla+0.1;
	$xy_file.="$stlo $stla\n";
	$t_file.="$stlo $stla 7 0 4 CB $stas[$k]\n";
      }

      if (defined $dfile and defined $sfile) {
	$data_files.="$dfile 0 $k\n"; $syn_files.="$sfile 0 $k\n";
	if ($plot_new) {$new_syn_files.="$nsfile 0 $k\n";}
	# windows
	if ($flexwin) {
	  @wins=split("\n",$file{$stas[$k]}{"win-$comps[$i]"});
	  for ($jj=0;$jj<@wins;$jj++) {
	    ($tstart,$tend)=split(" ",$wins[$jj]);
	    $tstart-=$file{$stas[$k]}{P}; $tend-=$file{$stas[$k]}{P}; #align t1
	    $k1=$k-0.5; $k2=$k+0.5;
	    $win_files.="$tstart $k1\n$tend $k1\n$tend $k2\n$tstart $k2\n>\n";}
	}
      }
    }
    chomp($data_files); chomp($syn_files); chomp($new_syn_files);

    # plot windows
    if ($flexwin) {
      chomp($win_files);chomp($win_files);
      plot_psxy(\*BASH,$psfile,"$JX $R $B[$i] -L -m $gwin ","$win_files");}

    # plot data/syn
    plot_pssac2_raw(\*BASH,$psfile,"$JX $R $B[$i] $E $C -M$size/0 $wdata -N","$data_files");
    plot_pssac2_raw(\*BASH,$psfile,"-JX -R -B $E $C -M$size/0 $wsyn -N","$syn_files");
    if ($plot_new) {
      plot_pssac2_raw(\*BASH,$psfile,"-JX -R -B $E $C -M$size/0 $wnsyn -N","$new_syn_files");
    }
    if ($i==0) {
      chomp($tex_file); chomp($t_file); chomp($xy_file);
      plot_pstext(\*BASH,$psfile,"-JX -R -B","$tex_file");}

    # add window properties
#    if (defined $mdir) {}

    plot_psxy(\*BASH,$psfile,"-JX -X-$x -Y-$y","");
  } # end of components

# add basin map here
  $rlat0=$file{$stlas[0]}{stla}-0.1; $rlat1=$file{$stlas[-1]}{stla}+0.1; 
  $rlon0=$file{$stlos[0]}{stlo}-0.1; $rlon1=$file{$stlos[-1]}{stlo}+0.1;
  $region="$rlon0/$rlon1/$rlat0/$rlat1";
  $BB=ceil(($rlon1-$rlon0)/10.);
  if ($socal) {
    plot_sc_all(\*BASH,$psfile,"$JM -W1 -w0.5 -R$region $jxy_map -B$BB/${BB}WesN","");}
  else{
    plot_pscoast(\*BASH,$psfile,"$JM  -R$region $jxy_map -B$BB/${BB}WesN");}
  if (defined $opt_m) {
    plot_psmeca(\*BASH,$psfile,"-JM -R -Sm0.2 ","$cmt_file");}
#  plot_psxy(\*BASH,$psfile,"-JM -R -St0.06 -W1 -G0/0/255","$xy_file");
  plot_pstext(\*BASH,$psfile,"-JM -R -B -G0/0/255","$t_file");

  plot_psxy(\*BASH,$psfile,"-JX -O","");
} # end of each plot

print BASH "ps2pdf.bash ds_??.ps\npdcat -r ds_??.pdf ds_all.pdf\n";

close(BASH);

system("chmod a+x plot_ds.bash; plot_ds.bash; rm -f ds_??.pdf ds_??.ps");
