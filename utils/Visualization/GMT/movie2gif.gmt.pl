#!/usr/bin/perl -w
#
# (version uses xcreate_movie_shakemap_AVS_DX_GMT)
#
# copy this script into your bin/ folder and run for example:
#
#  /movie2gif.gmt.pl -R -118.9/-117.9/33.1/34.1 -g -f 1/10 -n  -x
#
# will create gif files ../OUTPUT_FILES/gmt_movie_000***.gif

use Getopt::Std;

# USERS: point this to your SPECFEM3D/utils/lib directory
use lib '../../../utils/lib';
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
      to get clear/smoother movie images, try set longer hdur in the simulation.

EOF
exit(1);
}

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
  system("filelist.pl -I $sta_file -f /opt/seismo-util/data/STATIONS_TRINET  -n sta.tmp");}

if (not defined $opt_n) {$opt_n = 1;}

if ($opt_R) {$R = "-R$opt_R";} else {$R = "-R-120.3/-114.7/32.2/36.8";} #-120.3
$R2 = "${R}/0/0.8";
$R3 = "${R}/-1000/4000";
$JM = "-JM4.5i";
$JM3 = "${JM} -JZ1.5i";
$B = " -B2/2SWen";
$B2 = "-B2/2/0.5WSen";
$B3 = "-B2/2/1000WSen";

#$Itopo = "-I1m/1m";
$I = "-I0.5m";
if ($opt_d) {$S = "-S${opt_d}k";} else {$S = "-S5k";}

# USERS: adapt to your system
$datalib = "/opt/seismo-util/data/datalib_SC";
$fault_file = "$datalib/jennings.xy";
$fault_file_3d = "$datalib/jennings.xyz";
$cpt_file = "disp.cpt";
$topo_cpt = "$datalib/topo_cal_color.cpt";
$csh_file = "movie2gif.csh";
$big_grd_file = "$datalib/big_sc.grd";
$bin = ".";
$outdir = "../OUTPUT_FILES";

$factor = 10000;
$power = 1; #no distortion
$fac_max = 4;

# *******************
#   convert from movie data to gmt xyz files
if ($tran_to_gmt) {

  if (not -f "../DATA/Par_file") {die(" Check if ../DATA/Par_file exist or not\n");}
  (@junk) = split(" ", `grep '^DT' ../DATA/Par_file`);
  $dt = $junk[2];

  (@junk) = split(" ", `grep '^NTSTEP_BETWEEN_FRAMES' ../DATA/Par_file`);
  $nframe = $junk[2];

  # start and end time
  $startt = $start * $nframe;
  $endtt = $end * $nframe;

  print "Transfer the movie data files to gmt xyz files \n";
  system("$bin/xcreate_movie_shakemap_AVS_DX_GMT <<EOF > movie.log\n 3\n $startt\n $endtt\n 1\n 1\n EOF\n");
  if ($? != 0) {die("Error create movie xyz files\n"); }


} else {

  if (not -f "../DATA/Par_file") {die(" Check if ../DATA/Par_file exist or not\n");}
  (@junk) = split(" ", `grep '^DT' ../DATA/Par_file`);
  $dt = $junk[2];

  (@junk) = split(" ", `grep '^NTSTEP_BETWEEN_FRAMES' ../DATA/Par_file`);
  $nframe = $junk[2];

}
print "\ndt: $dt\nnframe: $nframe\n";


# figure out the maximum value
$max = `grep 'maximum amplitude' movie.log | awk 'BEGIN{max=0.0}{if(\$7 > max) max = \$7;}END{print max;}'`;
$max = $max / $fac_max;
print " The maximum value of all frames is : $max \n";

# ********************* write c-shell file
open(CSH,">$csh_file");

print CSH "gmtset PLOT_DEGREE_FORMAT D BASEMAP_TYPE  plain  ANNOT_FONT_SIZE_PRIMARY  12 HEADER_FONT_SIZE 15 \n";
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
  convert_linear("temp1.cpt",$cpt_file,$power);
} else {
  if (not -f $big_grd_file) {die("No such $big_grd_file \n");}
  print CSH "grdcut $big_grd_file $R -Gtopo.grd \n";
}

foreach $frame ($start .. $end) {
  if ($frame < 10) {$frame = sprintf("0%02d",$frame);}
  if ($frame < 100) {$frame = sprintf("%03d",$frame);}
  $file = "$outdir/gmt_movie_000$frame";
  $time = $frame * $nframe * $dt;
  $xyz_file = "$file.xyz"; if (not -f $xyz_file) {die("check $xyz_file\n");}
  $ps_file = "$file.ps";
  $gif_file = "$file.gif";
#  $tran_gif_file = "${file}_tran.gif";
  $grd_file = "$file.grd";

  print CSH "\n# frame $frame\n";

  if ($near) {
    # uses a nearest neighbor interpolation
    print CSH "nearneighbor -F $R $I $S -G$grd_file $xyz_file -E0 \n";
  }else{
    # surface gridding
    print CSH "surface $xyz_file $R $I -G$grd_file -T0.25 -C0.1 \n";
  }

#  else {
##    print CSH "awk '{print \$1,\$2,\$3/$max}' $xyz_file | pscontour $R $JM $B -A- -C$cpt_file -I -P> $ps_file\n";}  next;

  # divides grid by maximum value
  print CSH "grdmath $grd_file $max DIV = field1.grd\n";


  if ($no_topo) {

    # no topography

    print CSH "grdgradient field1.grd -V -A225 -Nt -Gfield.int\n";

    # creates image
    if ($two_d) {
      # grdimage for no topo
      # 2-d image
      print CSH "grdimage $JM $R field1.grd -Ifield.int  -C$cpt_file -K -P -V> $ps_file\n";
    } else {
      # grdview for no topo
      # 3-d image
      print CSH "psbasemap -P $E $R3 $B2 $JM3  -K > $ps_file\n";
      print CSH "grdview field1.grd -Qs -P $E $R2 $JM3 -C$cpt_file -Ifield.int -K -O >> $ps_file\n";
    }

  } else {

    # uses topography

    print CSH "grdmath field1.grd $factor MUL topo.grd ADD = field2.grd \n";
    print CSH "grdgradient field2.grd -V -A190/60 -Nt -Gfield.int\n";

    # creates image
    if ($two_d) {
      # for grdimage topo
      # 2-d image
      print CSH "grdimage field2.grd $JM $R -C$topo_cpt $B -Ifield.int -K -P  > $ps_file\n";}
    else {
      # grdview for topo
      # 3-d image
      print CSH "psbasemap -P $E $R3 $B3 $JM3 -K > $ps_file\n";
      print CSH "grdview field2.grd -Qs -P -Gtopo.grd $E $R3 $JM3 -C$topo_cpt -Ifield.int -K -O >> $ps_file \n";}
  }


  # coast lines, faults...
  if ($two_d) {

    # 2-d image add-ons

    print CSH "pscoast $JM $R -W4 -Na -Dh -K -O -P -V -S205/255/255 >> $ps_file \n";
    print CSH "psxy $JM $R -M -W2 -K -O -P -V $fault_file >> $ps_file\n";
    if ($opt_m) {print CSH "psxy $JM $R -Sa0.15 -W1 -G255/0/0 -K -O -P -V <<EOF >> $ps_file\n$elon $elat\nEOF\n";}
    if ($opt_s) {
      print CSH "awk '{print \$4, \$3}' sta.tmp | psxy $JM $R -St0.12 -W0.8 -G0/255/0 -K -O -P -V >> $ps_file\n";
      print CSH "awk '{print \$4, \$3+0.1, 12, 0, 4, \"CM\", \$1}' sta.tmp | pstext $JM $R -G0/0/255 -N -P -K -O -V >> $ps_file \n";
    }
    print CSH "pstext $JM $R -N -W255/255/255 -O  -P -V <<EOF >>$ps_file \n -114.8 36.5 12 0 0 RT time = $time s \nEOF\n";
  } else {

    # 3-d image add-ons

    # coast lines
    print CSH "pscoast -JM -JZ -R -W5 -Na -Dh -K -O -P -V $E -S205/255/255 >>  $ps_file \n";
    print CSH "psclip -JM -JZ -R -K -O -P -Z0 >> $ps_file <<EOF\n-120.0 32.2\n-120.0 36.4 \n-114.7 36.4\n-114.7 32.2\n-120.0 32.2\nEOF\n";

    # faults
    #print CSH "psxyz -JM -JZ -R -M -W2 -K -O -P -V $E $fault_file_3d >> $ps_file\n";
    print CSH "psclip -C -K -O >> $ps_file\n";

    # cmt solution
    if ($opt_m) {
      print CSH "psxyz -JM -JZ -R -Sa0.15 -W1 -G255/0/0 -K -O -P -V $E <<EOF >> $ps_file\n$elon $elat 0 \nEOF\n";
    }

    # stations names
    if ($opt_s) {
      print CSH "awk '{print \$4, \$3, 0}' sta.tmp | psxyz $E -JM -JZ -R -St0.12 -W0.8 -G0/255/0 -K -O -P -V >> $ps_file\n";
      print CSH "awk '{print \$4, \$3+0.15, 15, 0, 4, \"CM\", \$1}' sta.tmp | pstext -JM -JZ -G0/0/255 $E -R -N -P -K -O -V -Z0>> $ps_file \n";
    }

    # text for time
    print CSH "pstext -JM -JZ -R -N -W255/255/255 -O -P -V -Z0<<EOF >>$ps_file \n -114.8 36.5 12 0 0 RT time = $time s \nEOF\n";
  }

  # converts ps to gif
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
    }
    elsif($line =~ /^-/){ # fix bugs when data range is negative
      ($v1,$r1,$g1,$b1,$v2,$r2,$g2,$b2) = split(" ",$line);
      $v1 = -(abs($v1) ** $power), $v2 = -(abs($v2) ** $power);
      printf CPT2  ("%-15.10f\t%s\t%s\t%s\t%-15.10f\t%s\t%s\t%s\n",
                    $v1,$r1,$g1,$b1,$v2,$r2,$g2,$b2);
    }else {print CPT2 "$line";}
  }
  close(CPT2);
}
