#!/usr/bin/perl -w

use lib '/opt/seismo/lib/perl';
use GMT_PLOT;
use CMT_TOOLS;
use GMT_PLOT_SC;
use POSIX;

# specify the grid search parameter range
if (@ARGV == 0 or @ARGV < 2) {die(" plot_grid.pl cmt_file grid-output-files \n");}

$cmt=$ARGV[0];
if (not -f $cmt) {die("Check $cmt file\n");}
(undef,undef,$ename)=split(" ",`grep 'event name' $cmt`);
($sdr)=split(" ",`cmtsol2faultpar.pl $cmt | grep strike | awk '{print \$2}' | cut -d / -f 1-3`);
($sdr2)=split(" ",`cmtsol2faultpar.pl $cmt | awk 'NR==9 {print \$2}' | cut -d / -f 1-3`);
print "CMT file -- $ename ($sdr, $sdr2)\n";

# this controls the paper orientation;
$GMT_PLOT::paper_orient = "-P";

# global plotting pars
$nrows= 6;
$ncols= 3;
@name=("S","D","R");
for ($i=1;$i<=$nrows-1;$i++) {
  $Bs[$i]="Wesn"; $Bd[$i]="wesn"; $Br[$i]="wEsn";}
$Bs[0]="WesN"; $Bd[0]="wesN"; $Br[0]="wEsN";
$Bs[-1]="WeSn"; $Bd[-1]="weSn"; $Br[-1]="wESn";

# plot size jx/jy, plot origin in xy
($jx,$jy) = shift_xy("1 1.5","$ncols $nrows","7  8", \@xy,"0.84 0.84");
$JX = "-JX$jx/$jy";
$GMT_PLOT::paper_orient="-P";


# loop over input files
for ($k=1;$k<@ARGV;$k++) {

  $result=$ARGV[$k];
  print " == file = $result \n";
  if (not -f $result) {die("Check if $result exist or not\n");}

  $bash_file="plot_$k.bash";
  open(CSH,">$bash_file");
  print CSH "gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 9 HEADER_FONT_SIZE 10 MEASURE_UNIT inch PAPER_MEDIA letter TICK_LENGTH 0.1c\n";
  $psfile = "grid_$k.ps";

  ($ss,$es,$ns,$ts,$sd,$ed,$nd,$td,$sr,$er,$nr,$tr,$min0,$max0,@array) = pick_ds_dd_dr($result,$nrows);
  # LQY -- hard-wire this maximum color control
  if ($max0 > 100*$min0) {$max0=100*$min0;}

  $j=0;
  @as=@array[$j..$j+$ns-1]; $j+=$ns;
  @ad=@array[$j..$j+$nd-1]; $j+=$nd;
  @ar=@array[$j..$j+$nr-1]; $j+=$nr;
  @ps=@array[$j..$j+$nrows-1]; $j+=$nrows;
  @pd=@array[$j..$j+$nrows-1]; $j+=$nrows;
  @pr=@array[$j..$j+$nrows-1]; $j+=$nrows;

  print "s,d,r=$ss,$es,$ns,$ts,  $sd,$ed,$nd,$td,  $sr,$er,$nr,$tr\n";
  print "min/max=$min0/$max0\n";
  print "as=@as\nad=@ad\nar=@ar\n";
  print "ps=@ps\npd=@pd\npr=@pr\n";

 # plotting parameters
  @R = ("-R$sd/$ed/$sr/$er","-R$ss/$es/$sr/$er","-R$ss/$es/$sd/$ed");
  $min = $min0 - ($max0-$min0)/100;
  $max = $max0 + ($max0-$min0)/100;
  # plot intervals for psscale
  $db = floor(($max-$min)/4 * 1000)/1000;
  $db2 = $db/5; # makecpt -T interval
#  $db3 = $db/4; #grdcontour annotation interval -A
  $db4 = $db/4; # grdcontour interval -C
  $A = "-G1.2 -C$db4"; # grdcontour -A${db3}f10
  $db = "-B$db"; # psscale
  print "db = $db \n";
  # annotation intervals for psbasemap
  $dds = floor(($es-$ss)/5);
  $ddr=floor(($er-$sr)/5);
  $ddd=floor(($ed-$sd)/5);

  # plot position of the text for dip-rake, strike-rake, and strike-dip plot
  @pos= ("$ad[1] $ar[-2]","$as[1] $ar[-2]","$as[1] $ad[-2]");
  $arr=$ar[-1]+$ddr*2; # position for title
  $pp="$as[$nrows/2] $arr";

  print CSH "makecpt -T$min/$max/$db2 -Cseis -Z > temp.cpt \n";

  plot_psxy(\*CSH,$psfile,"$JX -R0/1/0/1 -K -X0 -Y0","");

  for ($i=0;$i<3;$i++) { #ij column
    for ($j=0;$j<$nrows;$j++) { #ii row
      $xy = $xy[$j][$i]; print "--col=$i, row=$j,  xy = $xy\n";
      ($x,$y) = split(/\//,$xy);
      plot_psxy(\*CSH,$psfile,"$JX -X$x -Y$y","");
      if ($i==0) {
  $value=$ps[$j]; $B="-B$ddd/${ddr}$Bs[$j]";
        if ($j==$nrows-1) {$B="-B$ddd:\"Dip\":/${ddr}:\"Rake\":$Bs[$j]";}
  print CSH "awk '\$1==$ps[$j] {print \$2,\$3,\$4}' $result > tmp1\n";
  print CSH "xyz2grd tmp1 -Gout.grd -I$td/$tr -R$sd/$ed/$sr/$er \n";
      } elsif ($i==1) {
  $value=$pd[$j]; $B="-B$dds/${ddr}$Bd[$j]";
        if ($j==$nrows-1) {$B="-B$dds:\"Strike\":/${ddr}$Bd[$j]";}
  print CSH "awk '\$2==$pd[$j] {print \$1,\$3,\$4}' $result > tmp1\n";
  print CSH "xyz2grd tmp1 -Gout.grd -I$ts/$tr -R$ss/$es/$sr/$er \n";
      } else {
  $value=$pr[$j];$B="-B$dds/${ddd}$Br[$j]";
        if ($j==$nrows-1) {$B="-B$dds:\"Strike\":/${ddd}:\"Dip\":$Br[$j]";}
  print CSH "awk '\$3==$pr[$j] {print \$1,\$2,\$4}' $result > tmp1\n";
  print CSH "xyz2grd tmp1 -Gout.grd -I$ts/$td -R$ss/$es/$sd/$ed \n";
      }
      print CSH "grdimage out.grd $R[$i] -JX $B -Ctemp.cpt -K -O -P >> $psfile \n";
      print CSH "grdcontour out.grd -R -JX -K -O -P $A>> $psfile\n";
      plot_pstext(\*CSH,$psfile,"-JX -W255 ","$pos[$i] 8 0 4 LT $name[$i]=$value");
      if ($i==1 and $j==0) {plot_pstext(\*CSH,$psfile,"-JX -W255 -N","$pp 14 0 4 LT $ename($sdr, $sdr2)");}
      plot_psxy(\*CSH,$psfile,"-JX -X-$x -Y-$y","");
      print CSH "\n";
    }
  }
print CSH "psscale -Ctemp.cpt -D4i/0.6i/2.5i/0.17h $db -K -O -V -P >> $psfile\n";
plot_psxy(\*CSH,$psfile,"$JX -O","");
close(CSH);

print "Executing $bash_file\n";
system("bash $bash_file");

}



sub pick_ds_dd_dr {

  my($file,$nrows) = @_;
  open(FILE,"$file") || die("Check if $file exists or not\n");
  @lines = <FILE>;
  close(FILE);
  %nns=(); %nnd=(); %nnr=(); @ps=(); @pd=(); @pr=();
  ($ss,$es,$sd,$ed,$sr,$er,$min,$max) = split(" ",`minmax -C $file`);

  for $line (@lines) {
    ($s,$d,$r,$v) = split(" ",$line);
    if ($v == $min) {$mins=$s;$mind=$d;$minr=$r;}
    $nns{$s}++; $nnd{$d}++; $nnr{$r}++;
  }
  $ns=keys %nns;
  $nd=keys %nnd;
  $nr=keys %nnr;
  @as=sort { $a <=> $b } keys %nns;
  @ad=sort { $a <=> $b } keys %nnd;
  @ar=sort { $a <=> $b } keys %nnr;
  $ts=$as[1]-$as[0]; $td=$ad[1]-$ad[0]; $tr=$ar[1]-$ar[0];
  if ($nrows > @as or $nrows > @ad or $nrows > @ar) {
      die("Choose a smaller nrows to be in accordance with ns,nd,nr\n");}

  $nleft=floor($nrows/2);
  $nright=$nrows-$nleft;
  $ns1=floor($ns/$nrows); $nd1=floor($nd/$nrows); $nr1=floor($nr/$nrows);
  $is=int(($mins-$ss)/$ts); $id=int(($mind-$sd)/$td); $ir=int(($minr-$sr)/$tr);

  for ($i=-$nleft;$i<$nright;$i++) {
    if ($is+$i*$ns1 >= @as) {$iss=$is+$i*$ns1-@as;} else {$iss=$is+$i*$ns1;}
    if ($id+$i*$nd1 >= @ad) {$idd=$id+$i*$nd1-@ad;} else {$idd=$id+$i*$nd1;}
    if ($ir+$i*$nr1 >= @ar) {$irr=$ir+$i*$nr1-@ar;} else {$irr=$ir+$i*$nr1;}

    @ps=(@ps,$as[$iss]);
    @pd=(@pd,$ad[$idd]);
    @pr=(@pr,$ar[$irr]);
  }

  return ($ss,$es,$ns,$ts,$sd,$ed,$nd,$td,$sr,$er,$nr,$tr,$min,$max,@as,@ad,@ar,@ps,@pd,@pr);


}
