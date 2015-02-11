#!/usr/bin/perl -w

#
#
#  If using CMT files, you must uncomment the opt_m commands below.
#
#  EXAMPLES (execute from main directory):
#  STATIONS file is produced in make_sac_files.pl
#
#  --> See multi_plot_rs.pl <--
#
#  scripts/plot_rs.pl -R0/200/0/200 -Ekt0 -M0.0003/s -W0.5p OUTPUT/*sacdat
#
#  scripts/plot_rs.pl -R0/200/0/200 -Ekt0 -m socal_adjoint_stf -n "Socal adjoint source time functions" -M0.0003/s -W0.5p OUTPUT/*sacstfadj
#
#  scripts/plot_rs.pl -R0/200/0/200 -Ekt0 -M0.0003/s -W0.5p -w0.5/255/0/0tap -S OUTPUT,sacsyn -D OUTPUT,sacdat,1,STATIONS
#

use Getopt::Std;
use lib '/opt/seismo-util/lib/perl';
use SACLST;
use SAC_TOOLS;
use CMT_TOOLS;
use GMT_PLOT_SC;
use GMT_PLOT;
use FILE_TOOLS;
use vars qw($opt_s $opt_r);

sub Usage(){

print STDERR <<EOF;
 Usage: plot_rs.pl

         this script plots record section for given data and synthetics,
         limit stations by azimuth (-A) and amplitude ratio (-T)

         -m : psfile label (formerly CMTSOLUTION)
         -n : title for the plot
         -R : default [-20/180/0/320]
         -B : tick marks, default [20/20]
         -C : cutting of trace, default no cutting
         -J : projection, default [-JX6/-8]
         -E : alignment of traces , like -Ekt5 default [-Ek8.0]
   -M : scale of trace, default [0.03/1]
         -W : data trace property
   -w : syn trace property
         -d data_file, data_ext, comp
         -D data_dir, data_ext, comp, station_file
         -S syn_dir, syn_ext, shift_hdr
         -A az_min/az_max
         -s plot LA basin instead of the whole SC map
         -r add text to the right instead of above the station in map.ps
         -T : maximum amplitude ratio of data/syn to plot
         data_files

EOF
exit(1);
}

@ARGV>0 or Usage();
if(!getopts('R:B:C:J:e:E:M:W:w:d:D:S:A:m:n:srT:')) {die("Check the input options\n");}
if ($opt_R) {$R = "-R$opt_R";} else {$R = "-R-20/180/-10/320";}
if ($opt_B) {$B =   "$opt_B";} else {$B = "20/20";}
if ($opt_C) {$C = "-C$opt_C";} else {$C = "";}
if ($opt_J) {$J = "-J$opt_J";} else {$J = "-JX6/-8";}
if ($opt_E) {$E = "-E$opt_E";} else {$E = "-Ek8.0";}
if ($opt_M) {$M = "-M$opt_M";} else {$M = "-M0.03/1";}
if ($opt_W) {$W = "-W$opt_W";} else {$W = "-W5";}
if ($opt_w) {$w = "-W$opt_w";} else {$w = "-W5/255/0/0";}
if ($opt_d) {
  ($data_file,$data_ext,$comp) = split(/\,/,$opt_d);
  $data_ext = $data_ext ? ".$data_ext" : "";
  $comp = $comp ? $comp : "BHT";
  @tmp=list2files($data_file);
  @data_files = ();
  foreach $data (@tmp) {
    (undef,$comp1) = split(" ",`saclst kcmpnm f $data`);
    if ($comp1 ne $comp) {next;}
    push(@data_files,"${data}${data_ext}");}}

elsif ($opt_D) {
  ($data_dir,$data_ext,$comp,$station_used) = split(/\,/,$opt_D);
  if (not $comp ) {die("need to specfiy component to plot\n");}
  if (not $data_dir) {die("need to specify the data directory\n");}
  if (not $data_ext) {die("need to specify the data extensions\n");}
  @data_files = ();
  if ($station_used) {
    (@sta) = list2files($station_used);
    foreach $sta (@sta) {
      @files = glob("$data_dir/*$sta.*$comp*$data_ext");
      #print "\n $sta $#files $files[0]";
      if (@files == 1) {push(@data_files,$files[0]);}
    }
  }else {
    @data_files = glob("$data_dir/*$comp*$data_ext");}
}else {@data_files = @ARGV; $comp = "BH";}

$numd = $#data_files + 1;
print "\n $numd data files, first is $data_files[0] \n";

if ($opt_A) {($min_az,$max_az) = split(/\//,$opt_A);
@data_files = sac_files_prune("az",$min_az,$max_az,@data_files);}

if (@data_files == 0) {die("No data files available to plot\n");}
@data_files = sac_files_sort("dist",1,@data_files);

if($opt_S) {
  $plot_syn=1; ($syn_dir,$syn_suff,$syn_shift_hdr) = split(/\,/,$opt_S);
  if (defined $syn_shift_hdr) {$S = "-Su$syn_shift_hdr"; }else {$S = "";}}
else {$plot_syn = 0;}
if($plot_syn & $M !~/\//){die "Don't use relative scale when plotting both data and synthetics -- change -M option \n";}

#------------------------

#if ($opt_m) {$cmt_file = $opt_m; ($evnm) = get_cmt_evnm($cmt_file);}
#else {$evnm = "socal_seis";}

if ($opt_m) {$evnm = $opt_m;}
else {$evnm = "socal_seis";}

if ($opt_n) {$title = $opt_n;}
else {$title = "Record section for SoCal";}

$B = "-B$B:.\"${title}\":WeSN";

#$psfile = "${evnm}.${comp}.ps";

# include the axis limits in the name of the figures
($junk,$R5) = split("-R",$R);
($R1,$R2,$R3,$R4) = split(/\//,$R5);
$pslab = $R1."_".$R2."_".$R3."_".$R4;

# shift y-limits by some fraction
$frac = 0.05;
$R3 = $R3 - $frac*($R4-$R3);
$R4 = $R4 + $frac*($R4-$R3);
$R  = "-R/$R1/$R2/$R3/$R4";

$name1 = ${evnm}."_".$pslab;
$psfile = $name1.".ps";
$jpgfile = $name1.".jpg";

print "\nFile labeling: \n";
print "evnm    = $evnm \n";
print "title   = $title \n";
print "psfile  = $psfile \n";
print "jpgfile = $jpgfile \n";

$name2 = "socal_map";
$psfile2 = $name2.".ps";
$jpgfile2 = $name2.".jpg";
$csh_file = "rs.csh";

#------------------------

open(CSH,">$csh_file") or die("Unable to open $csh_file\n");
print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch PLOT_DEGREE_FORMAT D BASEMAP_TYPE plain HEADER_FONT_SIZE 20\n";

plot_psxy(\*CSH,$psfile,"$J $R $B -P -K", "");

# plot the synthetics
# CAREFUL: LABELING SCHEME FOR DATA AND SYNTHETICS MUST FIT A REQUIRED FORMAT
if ($plot_syn) {
  ($temp1,$temp2) = match_data_and_syn("$syn_dir,$syn_suff","@data_files","");
  (@data_files) = split(" ",$temp1);
  (@syn_files)  = split(" ",$temp2);

  $numd = $#data_files + 1; $nums = $#syn_files + 1;
  print "\n loop: $numd data files";
  print "\n loop: $nums syn files \n\n";

  if ($opt_T) {
    for ($i = 0; $i < @data_files; $i ++ ) {
      (undef,$datamax) = split(" ",`saclst depmax f $data_files[$i]`);
      (undef,$synmax) = split(" ",`saclst depmax f $syn_files[$i]`);
      if ($datamax > $synmax * $opt_T or $synmax > $datamax * $opt_T) {
  print "skip $data_files[$i] and $syn_files[$i] -- too large amp. ratio\n";}
      else {push(@datatmp,$data_files[$i]); push(@syntmp,$syn_files[$i]);}
    } @data_files = @datatmp; @syn_files = @syntmp;
  }
  if (@data_files == 0 or @syn_files == 0) {die ("No data or syn files to plot\n");}
  plot_pssac2(\*CSH,$psfile,"$R $C $J $E $M $w $S -P ","@syn_files");
}

# plot the data
plot_pssac2(\*CSH,$psfile,"$R $C $J $E $M $W $B -P ","@data_files");

# make and plot a map of stations

$station_file       = "RS_STATIONS";   open(OUT, ">$station_file");
$station_file_basin = "RS_SOCAL";      open(OUT2,">$station_file_basin");
$text = "";

foreach $data (@data_files) {

   # use saclst to extract various quantities (evel, evdp not accessible)
  (undef,$sta,$dist,$az,$stla,$stlo,$net)=split(" ", `saclst kstnm dist az stla stlo knetwk f $data`);

  #$text .= sprintf("$R1 %6.1f 10 0 4 BC %s\n $R2 %6.1f 10 0 4 BC %6.1f    %6.1f \n",$dist,$sta,$dist,$dist,$az);

  $R2p = $R2 + 0.01*($R2-$R1);                          # x-coordinate for label at right

  if($dist >= $R3 && $dist <= $R4)
  {$text .= sprintf("$R2p $dist 8 0 1 LM %5s %12.1f km %12.1f \n",$sta,$dist,$az);}

  $text2 .= "$stlo $stla\n";                            # create text for psxy call below

  if ($opt_s) { if ($opt_r) {$stlo = $stlo + 0.02;} else {$stla = $stla + 0.03;}}
  else {if ($opt_r) {$stlo = $stlo + 0.04;} else {$stla = $stla + 0.06;}}

  $text3 .= "$stlo $stla 8 0 4 BC $sta\n";              # create text for pstext call below

  # make output files for checking (originally for specfem3d.f90)
  print OUT "$sta\n";                                   # print station label to file
  print OUT2 sprintf("%5s %5s %14.6f %14.6f %14.6f %14.6f \n",$sta,$net,$stlo,$stla,0,0);  # no elevation, burial
}
chomp($text); chomp($text2); chomp($text3);
close(OUT);

plot_pstext(\*CSH,$psfile,"-G255/0/0 -N -JX -P",$text);            # labels at the right of the seismograms (rec.sec.)

plot_psxy(\*CSH,$psfile,"-JX -P -O","");                           # wrap up (rec.sec.)

if ($opt_s) {$R = "-R-119.2/-117.2/33.5/34.7";} else {$R = "";}

plot_sc_all(\*CSH,$psfile2,"-K $R -W4 -w3 ", "");

plot_psxy(\*CSH,$psfile2,"-St0.10 -G0/0/255 -W2",$text2);          # symbols at the stations (map)

plot_pstext(\*CSH,$psfile2,"-G255/0/0", $text3);                   # labels at the stations (map)

#if ($opt_m) {plot_psmeca(\*CSH,$psfile2,"-Sm0.2 -O",$cmt_file);}

#------------------------------------------------------------------------
close(CSH);
system("csh -fv $csh_file");
system("convert $psfile $jpgfile");
system("convert $psfile2 $jpgfile2");

#system("xv $jpgfile2 &"); # map with stations
#system("xv $jpgfile &");  # record section

#-----------------------------------------------------------------------

