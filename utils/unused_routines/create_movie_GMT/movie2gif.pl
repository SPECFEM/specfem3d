#!/usr/bin/perl -w

use Getopt::Std;
use lib '/opt/seismo-util/lib/perl';
use CMT_TOOLS;
use GMT_PLOT;
use vars qw($opt_t $opt_p $opt_n $opt_g $opt_x $opt_2 $opt_s $opt_d $opt_R);
sub Usage {
  print STDERR <<EOF;

      This script converts moviedata to xyz files, and then plot them in GMT
      Usage : movie2gif.pl -R -m cmt_file -g -f start/end -p -2 -n -d dist(km) -s sta_file -x
      -m cmt file
      -R min_lon/max_lon/min_lat/max_lat (default -120.3/-114.7/32.2/36.8)
      -g convert the movie data to gmt xyz files, you can choose to skip
         this step if you have already done it.
      -f start and end frame
      -p add topography to the plot
      -2 2d plot, otherwise 3d plot
      -s station name file to plot
      -n -- run nearneighbor command to interpolate xyz file into grd file, 
            since this is the step that takes most of the time, you can choose
            to skip this step if you have already run it.
      -d -- distance to average for nearneighbor command option -S (default 5 km)
      -x -- executing the c shell script

EOF
exit(1);
}	
#$bindir = "/opt/seismo-util/source/basin_inversion/bin";


if (@ARGV == 0) {Usage();}
if (not getopts('m:f:pg2xns:R:d:')) {die("Check options\n");}
if ($opt_m) {$cmt_file = $opt_m;
	     ($elat,$elon) = get_cmt_location($cmt_file);}

if (not defined $opt_f) {die("give the start and end of frame\n");}
($start,$end) = split(/\//,$opt_f);
#$transparent = 1;
if ($opt_p) {$no_topo = 0;} else {$no_topo = 1;}
if ($opt_g) {$tran_to_gmt = 1;} else {$tran_to_gmt = 0;}
if ($opt_2) {$two_d = 1;  $E = "";} else {$two_d = 0; $E = "-E190/60";}
if ($opt_n) {$near = 1;} else {$near = 0;}
if ($opt_s) {
  $sta_file = $opt_s;
  system("filelist.pl -I $sta_file -f $sta_file  -n sta.tmp");}

if (not defined $opt_n) {$opt_n = 1;}

if ($opt_R) {$R = "-R$opt_R";} else {$R = "-R-120.3/-114.7/32.2/36.8";} #-120.3
$R2 = "${R}/0/0.8";
$R3 = "${R}/-1000/4000";
$JM = "-JM4.5i";
$JM3 = "${JM} -JZ1.5i";
$B = " -B2/2wesn ";
$B2 = "-B2/2/0.5WSen";
$B3 = "-B2/2/1000WSen";

#$Itopo = "-I1m/1m";
$I = "-I0.5m";
if ($opt_d) {$S = "-S${opt_d}k";} else {$S = "-S5k";}

$datalib = "/opt/seismo-util/data/datalib_SC";
$fault_file = "$datalib/jennings.xy";
$fault_file_3d = "$datalib/jennings.xyz";
$cpt_file = "disp.cpt";
$topo_cpt = "$datalib/topo_cal_color.cpt";
$csh_file = "movie2gif.csh";
$big_grd_file = "$datalib/big_sc.grd";

$factor = 10000;


# *******************
#   convert from movie data to gmt xyz files
if ($tran_to_gmt) {
  print "Transfer the movie data files to gmt xyz files \n";
  $cmdline="xcreate_movie_GMT . Par_file $start $end . > movie.log";
  print "running command:\n   $cmdline ... \n";
  system($cmdline);
  if ($? != 0) {die("Error create movie xyz files\n"); }
  (@junk) = split(" ", `grep '^ DT' movie.log`); $dt = $junk[2];
  (@junk) = split(" ", `grep '^NTSTEP_BETWEEN_FRAMES' Par_file`); $nframe = $junk[2];
} else {
  if (not -f "movie.log") {die(" Check if movie log file exist or not\n");}
  $dt = 0.009; $nframe = 200;
}
# figure out the maximum value
(@junk) = split(" ", `grep 'maximum absolute value' movie.log`);
$max = $junk[-1] / 1;
print " The maximum value of all frames are : $max \n";

# *********************
open(CSH,">$csh_file");

print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 18 HEADER_FONT_SIZE 20\n";
# make cpt file, manipulate the white color part
if ($no_topo) {
  system("makecpt -T-1/1/0.1 -Cpolar  > temp.cpt");
  open(SED,">sed.in");
  print SED "/^N/s/^.*\$/N\t255\t255\t255/ \n";
  print SED "/^F/s/^.*\$/F\t255\t0\t0/ \n";
  print SED "/^B/s/^.*\$/B\t0\t0\t255/ \n";
  print SED "/^-0.1/s/^.*\$/-0.1\t255\t255\t255\t0\t255\t255\t255/ \n";
  print SED "/^5.55112e-17/s/^.*\$/0\t255\t255\t255\t0.1\t255\t255\t255/ \n";
  close(SED);
  system("sed -f sed.in temp.cpt > temp1.cpt\n");
  convert_linear("temp1.cpt",$cpt_file,3);
} else {
  if (not -f $big_grd_file) {die("No such $big_grd_file \n");}
  print CSH "grdcut $big_grd_file $R -Gtopo.grd \n";
}

foreach $frame ($start .. $end) {
  if ($frame < 10) {$frame = sprintf("00%02d",$frame);}
  if ($frame < 100) {$frame = sprintf("0%03d",$frame);}
  if ($frame < 1000) {$frame = sprintf("%04d",$frame);}
  $file = "gmt_movie_00$frame";
  $time = $frame * $nframe * $dt;
  $xyz_file = "$file.xyz"; if (not -f $xyz_file) {die("check $xyz_file\n");}
  $ps_file = "$file.ps";
  $gif_file = "$file.gif";
#  $tran_gif_file = "${file}_tran.gif";
  $grd_file = "$file.grd";

  print CSH "\n# frame $frame\n";
  if ($near) {
    print CSH "nearneighbor -F $R $I $S -G$grd_file $xyz_file -E0 \n";
    print CSH "grdmath $grd_file $max DIV = field1.grd\n";}
#  else {
##    print CSH "awk '{print \$1,\$2,\$3/$max}' $xyz_file | pscontour $R $JM $B -A- -C$cpt_file -I -P> $ps_file\n";}  next;

  if ($no_topo) {
    print CSH "grdgradient field1.grd -V -A225 -Nt -Gfield.int\n";
    if ($two_d) { # grdimage for no topo
      print CSH "grdimage $JM $R field1.grd -Ifield.int  -C$cpt_file -K -P -V> $ps_file\n";
    } else { # grdview for no topo
      print CSH "psbasemap -P $E $R3 $B2 $JM3  -K > $ps_file\n";
      print CSH "grdview field1.grd -Qs -P $E $R2 $JM3 -C$cpt_file -Ifield.int -K -O >> $ps_file\n"; }
  } else {
    print CSH "grdmath field1.grd $factor MUL topo.grd ADD = field2.grd \n";
    print CSH "grdgradient field2.grd -V -A190/60 -Nt -Gfield.int\n";
    if ($two_d) { # for grdimage topo
      print CSH "grdimage field2.grd $JM $R -C$topo_cpt $B -Ifield.int -K -P  > $ps_file\n";}
    else { # grdview for topo
      print CSH "psbasemap -P $E $R3 $B3 $JM3 -K > $ps_file\n";
      print CSH "grdview field2.grd -Qs -P -Gtopo.grd $E $R3 $JM3 -C$topo_cpt -Ifield.int -K -O >> $ps_file \n";}
  }
#  if ($no_topo) {
    if ($two_d) {
      print CSH "pscoast $JM $R -W4 -Na -Dh -K -O -P -V -S205/255/255 >> $ps_file \n";
      print CSH "psxy $JM $R -M -W2 -K -O -P -V $fault_file >> $ps_file\n";
      if ($opt_m) {print CSH "psxy $JM $R -Sa0.15 -W1 -G255/0/0 -K -O -P -V <<EOF >> $ps_file\n$elon $elat\nEOF\n";}
      if ($opt_s) {
	print CSH "awk '{print \$4, \$3}' sta.tmp | psxy $JM $R -St0.12 -W0.8 -G0/255/0 -K -O -P -V >> $ps_file\n";
	print CSH "awk '{print \$4, \$3+0.1, 12, 0, 4, \"CM\", \$1}' sta.tmp | pstext $JM $R -G0/0/255 -N -P -K -O -V >> $ps_file \n";
      }
    print CSH "pstext $JM $R -N -W255/255/255 -O  -P -V <<EOF >>$ps_file \n -114.8 36.5 12 0 0 RT time = $time s \nEOF\n";
  } else {
    print CSH "pscoast -JM -JZ -R -W5 -Na -Dh -K -O -P -V $E -S205/255/255 >>  $ps_file \n";
    print CSH "psclip -JM -JZ -R -K -O -P -Z0 >> $ps_file <<EOF\n-120.0 32.2\n-120.0 36.4 \n-114.7 36.4\n-114.7 32.2\n-120.0 32.2\nEOF\n";
    print CSH "psxyz -JM -JZ -R -M -W2 -K -O -P -V $E $fault_file_3d >> $ps_file\n";
    print CSH "psclip -C -K -O >> $ps_file\n";
    if ($opt_m) {print CSH "psxyz -JM -JZ -R -Sa0.15 -W1 -G255/0/0 -K -O -P -V $E <<EOF >> $ps_file\n$elon $elat 0 \nEOF\n";}
    if ($opt_s) {
      print CSH "awk '{print \$4, \$3, 0}' sta.tmp | psxyz $E -JM -JZ -R -St0.12 -W0.8 -G0/255/0 -K -O -P -V >> $ps_file\n";
      print CSH "awk '{print \$4, \$3+0.15, 15, 0, 4, \"CM\", \$1}' sta.tmp | pstext -JM -JZ -G0/0/255 $E -R -N -P -K -O -V -Z0>> $ps_file \n";
    }
    print CSH "pstext -JM -JZ -R -N -W255/255/255 -O -P -V -Z0<<EOF >>$ps_file \n -114.8 36.5 12 0 0 RT time = $time s \nEOF\n";
   }
#}
    print CSH "convert $ps_file $gif_file\n";
#  if ($no_topo and  $transparent) {print CSH "giftrans -t\\\#fefefe $gif_file > $tran_gif_file\n";}
}

close(CSH);
if ($opt_x) {system("csh -fv $csh_file");}


#---------------------------

sub convert_linear {

  my ($cpt,$new_cpt,$power) = @_;

  open(CPT,"$cpt"); @old_cpt = <CPT>;close(CPT);
  open(CPT2,">$new_cpt");
  foreach $line (@old_cpt) {
    if ($line =~/^\d+/) {
      ($v1,$r1,$g1,$b1,$v2,$r2,$g2,$b2) = split(" ",$line);
      $v1 = $v1 ** $power, $v2 = $v2 ** $power;
      printf CPT2  ("%-15.10f\t%s\t%s\t%s\t%-15.10f\t%s\t%s\t%s\n",
                    $v1,$r1,$g1,$b1,$v2,$r2,$g2,$b2);
    }else {print CPT2 "$line";}
  }
  close(CPT2);
}
