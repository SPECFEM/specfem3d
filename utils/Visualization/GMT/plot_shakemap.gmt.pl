#!/usr/bin/perl -w
#
# (version uses xcreate_movie_shakemap_AVS_DX_GMT)
#
# copy this script into your bin/ folder and run for example:
#
#  ./plot_shakemap.gmt.pl ../DATA/ 2 ../DATA/CMTSOLUTION
#
# will create shakemap file ../OUTPUT_FILES/gmt_shaking_map.ps for peak ground velocity

# USERS: change this to point to SPECFEM3D/utils/lib/
use lib '../utils/lib';
use GMT_PLOT;
use GMT_PLOT_SC;
use CMT_TOOLS;
use POSIX;

# usage
if (@ARGV != 3) {die("Usage: plot_shakemap.pl dir_name type CMTSOLUTION \n");}

$dir = $ARGV[0]; $type = $ARGV[1]; $cmt_file = $ARGV[2];
if (not -d $dir or not -f $cmt_file) {die("not file $dir or $cmt_file\n");}
if ($type != 1 and $type != 2 and $type != 3) {die("type can only be 1,2,3\n");}
@name = ("Displacement", "Velocity", "Acceleration");

$power = 0.5;
$file = "../OUTPUT_FILES/gmt_shaking_map";
$xyz_file = "${file}.xyz";
$grd_file = "${file}.grd";
$cpt_file = "${file}.cpt";

# USERS: take your own input grid file:
$int_file = "/opt/seismo-util/data/datalib_SC/topo_cal.int";

$colorbar = "rainbow";
$bindir = ".";
$title .= "Peak $name[$type-1]";

if (not -f "DATA/Par_file") {
  if (not -f "$dir/Par_file") {die("Check if $dir/Par_file exist or not\n");}
  system("mkdir -p DATA/; cp $dir/Par_file DATA/Par_file");}

($elat,$elon) = get_cmt_location($cmt_file);

print "Extract the shaking data ...\n";
system("$bindir/xcreate_movie_shakemap_AVS_DX_GMT <<EOF > shake.out\n 3\n -1\n $type\n 1\n EOF\n");
if ($? != 0) {die("Error executing conversion to GMT files\n");}

# determines maximum value
(@tmp) = split(" ",`minmax -C $xyz_file`);
$max = $tmp[5];
print "\nmaximum: $max\n";

@units = ("m", "m/s","m/s^2");
$title.= sprintf(", Normalized by %5.2e $units[$type-1]",$max);

#normalizes data
`awk '{print \$1,\$2,\$3/max}' max=$max ${xyz_file} > ${xyz_file}.tmp ` ;
`mv ${xyz_file}.tmp ${xyz_file}`;

print "The output xyz file has been scaled to [0,1]\n";
$min = 0; $max = 0.5;
$nm = 50;
$inm = ($max-$min)/$nm;
system("makecpt -C$colorbar -T$min/$max/$inm -V -Z> $cpt_file\n");
convert_linear($cpt_file,"$cpt_file.1",1/$power);


$GMT_PLOT::paper_orient = "-P";
$R = "-118.9/-117.9/33.1/34.1";$I1 = 0.00833333;
($lon1,$lon2,$lat1,$lat2) = split(/\//,$R);
$nlon = int(($lon2-$lon1)/$I1 +0.5); $nlat = int(($lat2-$lat1)/$I1 + 0.5);
$R = "-R$R";

$I = "-I1m/1m";
$S = "-S5k";
$J = "-JM5i";
$B = "-B1/1:.\"$title\":";

$csh_file = "plot_shakemap.csh";
$ps_file = "$file.ps";
open(CSH,">$csh_file");
print CSH "gmtset DEGREE_FORMAT = 3 BASEMAP_TYPE  = plain  ANNOT_FONT_SIZE_PRIMARY 12 HEADER_FONT_SIZE = 15\n";
print CSH "nearneighbor $R $I $S -G$grd_file -V ${xyz_file}\n";
print CSH "grdsample $R $grd_file -G$grd_file.1 -N$nlon/$nlat -F -V\n";
print CSH "grdcut $int_file -Gtopo_int.1 $R -V\n";

plot_psbasemap(\*CSH,$ps_file,"$R $J $B -K");
plot_grdimage(\*CSH,$ps_file,"-Itopo_int.1","$grd_file.1",$cpt_file);
plot_pscoast(\*CSH,$ps_file,"-Na -Dh");
plot_sc_faults(\*CSH,$ps_file,"-W4/200");
plot_psxy(\*CSH,$ps_file,"-Sa0.2","$elon $elat");
print CSH "psscale -D2.4i/-0.4i/3i/0.15ih -B0.05 -C$cpt_file.1 -P -K -O -V >>$ps_file\n";
plot_psxy(\*CSH,$ps_file,"-O","");
close(CSH);

system("csh -fv $csh_file");
print "DONE!!\n";


sub convert_linear {

  my ($cpt,$new_cpt,$power) = @_;

  open(CPT,"$cpt"); @old_cpt = <CPT>;close(CPT);
  open(CPT2,">$new_cpt");
  foreach $line (@old_cpt) {
    if ($line =~/^\d+/) {
      ($v1,$r1,$g1,$b1,$v2,$r2,$g2,$b2) = split(" ",$line);
      $v1 = $v1 ** $power, $v2 = $v2 ** $power;
      printf CPT2  ("%-10.8f\t%s\t%s\t%s\t%-10.8f\t%s\t%s\t%s\n",
		    $v1,$r1,$g1,$b1,$v2,$r2,$g2,$b2);
    }else {print CPT2 "$line";}
  }
  close(CPT2);
}
