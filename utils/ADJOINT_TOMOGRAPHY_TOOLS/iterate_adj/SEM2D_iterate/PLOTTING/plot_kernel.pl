#!/usr/bin/perl -w

#==========================================================
#
#  plot_kernel.pl
#  Carl Tape
#  15-Feb-2006
#
#  This script plots a kernel in GMT.
#
#  Example (no longer there):
#    scripts/plot_kernel.pl OUTPUT/run_1200/ OUTPUT/run_1200/event_005/ 25 132 3500 20.0 sr.txt fun_test0 -7/4
#    scripts/plot_kernel.pl OUTPUT/run_1200/ OUTPUT/run_1200/event_005/ 25 132 3500 20.0 sr.txt fun_test1 1/1
#    scripts/plot_kernel.pl OUTPUT/run_1200/ OUTPUT/run_1200/event_005/ 25 132 3500 20.0 sr.txt fun_test2 1/1
#    scripts/plot_kernel.pl OUTPUT/run_1200/ OUTPUT/run_1200/event_005/ 25 132 3500 20.0 sr.txt fun_test3 1/1
#
#    scripts/plot_kernel.pl OUTPUT 4480 5 16 1 132 3500 20.0 sr.txt fun 1/1
#==========================================================

if (@ARGV < 8) {die("Usage: plot_kernel.pl basedir iopt xxx \n");}
($dir1,$dir,$ievent,$iter,$nsrc,$nrec,$c0,$per,$srcrec,$name,$color) = @ARGV;
$comp = 1;

# color for the kernel
$scale_color = 21.0;
$colorbar = "seis";
@scale  = split("/",$color);
$pwr = $scale[0];
$norm = "1e${pwr}";
$cmax = "$scale[1]";
$ss = $cmax;
$ds = 2*$ss/$scale_color;
$bs = sprintf("%3.3e",0.9*$ss);  # colorbar
$TsK = sprintf("-T%3.3e/%3.3e/%3.3e",-$ss,$ss,$ds);
print "$color $norm $cmax \n TsK = $TsK \n";
$BscaleK  = sprintf("-B%2.2e:\" K\@-R\@- ( x, y, t )  ( 10\@+-%2.2i\@+  m\@+-2\@+ [\@~d\143\@~] )\":",$bs,$pwr);

#die("testing");

# file names
$mfile_syn = "socal_vel_syn.dat";
$rec_file  = "recs_lonlat.dat";
$src_file  = "events_lonlat.dat";

$file_msyn = "$dir1/$mfile_syn";
$file_rec  = "$dir1/$rec_file";
$file_src  = "$dir1/$src_file";
$file_src_rec = "$dir/$srcrec";
$file_ker  = "$dir1/$name.dat";

if (not -f $file_rec) { die("Check if $file_rec exist or not\n") }
if (not -f $file_src) { die("Check if $file_src exist or not\n") }
if (not -f $file_src_rec) { die("Check if $file_src_rec exist or not\n") }
if (not -f $file_ker) { die("Check if $file_ker exist or not\n") }
if (not -f $file_msyn) { die("Check if $file_msyn exist or not\n") }

$plabel = "/home/carltape/sem2d/2d_adjoint/scripts/plot_kernel.pl";

#die("\ntesting\n");

# plotting specifications
$fsize0 = "14";
$fsize1 = "11";
$fsize2 = "9";
$fsize3 = "5";
$fontno = "1";    # 1 or 4
$tick   = "0.2c";

# plot symbols for sources, receivers, and shelf
$src    = "-W1.0p/0/0/0 -Sa0.25";
$rec    = "-Sc10p -W0.5p";
$src0   = "-Sa0.25 -W1.0p/0/0/0 -G255";
$rec0   = "-Sc10p  -W1.0p/0/0/0 -G255";
#$rec    = "-W0.5p/0/0/0 -St0.10";
$Wshelf = "-W1.0/0/0/0tap";

#-------------------------
# color

# KEY: scaling for color
$scale_color = 21.0;
$colorbar = "seis";

#-------------------------

# write plotting scripts
$wid = 6;
$J1 = "-JM${wid}i";      # in lat-lon
$origin = "-X1 -Y3.5";

$B1syn = "-Ba1f0.25:.\" \":WESN";

$Dlen = 0.5*$wid;
$Dx = $wid/2;
$Dy = -0.40;
$Dscale1 = "-D$Dx/$Dy/$Dlen/0.20h";

#$bs1 = 0.5; $Bscale1  = sprintf("-B%2.2e:\" Phase Velocity ( km s\@+-1\@+ )\": -E10p",$bs1);
#$bs1 = 2; $Bscale1  = sprintf("-B%2.2e:\" Percent Perturbation from %3.3f  km s\@+-1\@+ \": -E10p",$bs1,$c0/1000);

#-------------------------

#$per = 20;
#$c0 = 3.5;
$per_lab = sprintf("T = %3.3f s",$per);
$c0_lab  = sprintf("c0 = %3.3f km/s",$c0/1000);
print "\n$per_lab, $c0_lab\n";

$title = "Rayleigh Wave Phase Velocity";

# set bounds for the plotting
#$xmin = -121; $xmax = -115.0; $zmin = 31.75; $zmax = 36.5;
($xmin,$xmax,$zmin,$zmax,$smin,$smax,$tmin,$tmax) = split(" ",`minmax -C $file_msyn`);
($tt) = sort{$b <=> $a} ($tt,abs($tmin),abs($tmax));
$dinc = 0.25;  # buffer around min/max
$xmin = $xmin-$dinc;  $xmax = $xmax+$dinc;
$zmin = $zmin-$dinc;  $zmax = $zmax+$dinc;
$R = "-R$xmin/$xmax/$zmin/$zmax";

print "\n $R \n $smin $smax \n -$tt $tt \n";

# plot title
$xtx = $xmin+0.5*($xmax-$xmin);
$ztx = $zmin+1.10*($zmax-$zmin);

#$name    = "kernel";
$psfile  = "$name.ps";
$jpgfile = "$name.jpg";

  #===============================================
  print "\nWriting CSH file...\n";
  $cshfile = "plot_kernel.csh";
  open(CSH,">$cshfile");
  print CSH "gmtset BASEMAP_TYPE plain TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 PLOT_DEGREE_FORMAT D HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1\n";
  #===============================================

  # make colorpoint file
  $cptfile1 = "color0.cpt";
  print CSH "makecpt -C$colorbar $TsK > temp1\n";
  print CSH "sed 's/^B.*/B       170     0      0  /' temp1  >  temp2\n";
  print CSH "sed 's/^F.*/F         0     0    255  /' temp2 > $cptfile1\n";

  $plot_title = "Model for synthetics";

  print CSH "psbasemap $B1syn $R $J1 -P -K -V $origin > $psfile\n";  # START

  # phase velocity map
  #print CSH "awk '{print \$1,\$2,\$3/1000}' $file_msyn | pscontour $R $J1 -A- -C$cptfile1 -I -P -K -O -V >> $psfile\n";
  print CSH "awk '{print \$1,\$2,\$3 / $norm}' $file_ker | pscontour $R $J1 -A- -C$cptfile1 -I -P -K -O -V >> $psfile\n";

  print CSH "pscoast $J1 $R $B1syn -W1p -Na/1p -Dh -P -K -O -V >> $psfile\n";
  #print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J1 $R $Wshelf -K -O -P -V >> $psfile\n";
  print CSH "psscale -C$cptfile1 $Dscale1 $BscaleK -K -P -O -V >> $psfile \n";

  # plot receivers with numbered label
  #print CSH "awk '{print \$1,\$2}' $file_rec |psxy -N $J1 $R -K -O -P -V $rec >> $psfile\n";
  print CSH "awk '{print \$1,\$2}' $file_rec | tail -n $nrec > tempk\n";
  print CSH "psxy tempk -N $J1 $R -K -O -P -V $rec >> $psfile\n";
  $file_rec2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_rec | tail -n $nrec > $file_rec2\n";
  print CSH "pstext $file_rec2 -N $J1 $R -K -O -P -V >> $psfile\n";

  # plot sources with numbered label
  print CSH "awk '{print \$1,\$2}' $file_src | tail -n $nsrc > tempv\n";
  print CSH "psxy tempv -N $J1 $R -K -O -P -V $src >> $psfile\n";
  #print CSH "awk '{print \$1,\$2}' $file_src |psxy -N $J1 $R -K -O -P -V $src >> $psfile\n";
  $file_src2 = text_src; $angle = 0; $just = "CM";
  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_src | tail -n $nsrc > $file_src2\n";
  print CSH "pstext $file_src2 -N -G255/0/0 $J1 $R -K -O -P -V >> $psfile\n";

#  # plot receivers with numbered label
#  print CSH "awk '{print \$1,\$2}' $file_rec |psxy -N $J1 $R -K -O -P -V $rec >> $psfile\n";
#  $file_rec2 = text_rec; $angle = 0; $just = "CM";
#  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_rec > $file_rec2\n";
#  print CSH "pstext $file_rec2 -N $J1 $R -K -O -P -V >> $psfile\n";

#  # plot sources with numbered label
#  print CSH "awk '{print \$1,\$2}' $file_src |psxy -N $J1 $R -K -O -P -V $src >> $psfile\n";
#  $file_src2 = text_src; $angle = 0; $just = "CM";
#  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_src > $file_src2\n";
#  print CSH "pstext $file_src2 -N -G255/0/0 $J1 $R -K -O -P -V >> $psfile\n";

  # plot title and GMT header
  $Utag = "-U/0/0.25/$plabel";
  $Ushift = "-Xa-$dX -Ya6.5";

  $plot_title = "Kernel";
  if (not -f $file_src_rec) { die("Check if $file_src_rec exist or not\n") }

  # plot individual receiver
  print CSH "awk '\$1 == \"R\" {print \$2,\$3}' $file_src_rec |psxy -N $J1 $R -K -O -P -V $rec0 >> $psfile\n";
  $file_rec2 = text_rec; $angle = 0; $just = "CM";
  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_rec | tail -n $nrec > $file_rec2\n";
  #print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_rec > $file_rec2\n";
  print CSH "pstext $file_rec2 -N $J1 $R -K -O -P -V >> $psfile\n";

  # plot individual source
  print CSH "awk '\$1 == \"S\" {print \$2,\$3}' $file_src_rec |psxy -N $J1 $R -K -O -P -V $src0 >> $psfile\n";
  $file_src2 = text_src; $angle = 0; $just = "CM";
  print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_src | tail -n $nsrc > $file_src2\n";
  #print CSH "awk '{print \$1,\$2,$fsize3,$angle,$fontno,\"$just\",\$3}' $file_src > $file_src2\n";
  print CSH "pstext $file_src2 -N $J1 $R -K -O -P -V >> $psfile\n";

#-----------------------------
  print CSH "pstext -N $J1 $R $Utag -O -P -V $Ushift >>$psfile<<EOF\n $xmin $zmin $fsize0 0 $fontno LM $plot_title\nEOF\n";  # FINISH

  print CSH "convert $psfile $jpgfile\n";

  close (CSH);
  system("csh -f $cshfile");
  system("xv $jpgfile &")

#=================================================================
