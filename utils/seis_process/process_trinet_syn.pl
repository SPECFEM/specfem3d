#!/usr/bin/perl


use Time::Local;
use Math::Trig;
use Getopt::Std;
use POSIX;
use lib '/opt/seismo-util/lib/perl/';
use CMT_TOOLS;


sub Usage{
print STDERR <<END;

Usage:   process_trine_syn_new.pl
     -S -m CMTFILE -h -o offset -lStart/End -tTmin/Tmax
     -i Dir -p -P n -a STAFILE -c  -d OutDir -x Ext
     -s sps -y t0/slope0/t1/slope1... name of files
where
    -S input files are already in sac format
    -m CMTfile name for event origin time,location and half duration
    -h convolve with triangle stf with half duration from -m in SAC
    -o moves the synthetics back by 'offset' from centroid time.
    -l lmin/lmax   start and end of trace
    -t Tmin/Tmax -- the shortest peroid of bandpass filter
    -i Dir -- convolve with instrument response in Dir (none for defaulted
      /opt/seismo-util/data/TriNet_PZs/)
    -p -- pick P and S arrivals according to socal model
    -P -- number of poles in butterworth filter (default 4)
    -a STAFILE -- add station information from STAFILE, use "-a none"
       to get information from Jascha''s station file STATIONS_TRINET
    -s sps -- resampling rate for synthetics ( default is 1.0)
    -y t0/slope0/t1/slope1 -- add to t3... the arrival time of waves
    -d OutDir -- directory to dump output (default: same as input files)
    -x Ext -- add extension Ext
    -c -- show sac output
    -w -- add CMTSOLUTION Moment information into header resp0 - resp5
          and the exponent information into isynth

    names of files -- syn files to be processed

END
exit(1);
}

@ARGV > 1 or Usage();

if (!getopts('Sm:ho:l:t:i:pP:a:cd:x:vs:y:w')) {die(" check input arguments\n");}
if ($opt_t) {($tmin, $tmax) = split(/\//,$opt_t);
             $f1 = 1./$tmax;  $f2=1./$tmin;}
if ($opt_x) {$ext=".".$opt_x;} else {$ext="";}
if (!$opt_s) {$dt=0.05;} else {$dt = 1.0/$opt_s;}
if (!$opt_P) {$opt_P=4;}
if ($opt_l) {($lmin,$lmax) = split(/\//,$opt_l);} else {$lmin = 0; $lmax = 3600;}
if ($opt_i and !-d $opt_i) {$opt_i="/opt/seismo-util/data/TriNet_PZs";}
if ($opt_a and !-f $opt_a) {$opt_a="/opt/seismo-util/data/STATIONS_TRINET";}
if ($opt_o and not $opt_m) {die("Specify centroid time first\n");}
if ($opt_d and not -d $opt_d) {die("No such directory as $opt_d\n");}


$saclst="/opt/seismo-util/bin/saclst";
$sacch="/opt/seismo-util/bin/sacch";
$phtimes="/opt/seismo-util/bin/phtimes";
$asc2sac="/opt/seismo-util/bin/ascii2sac.csh";
if (! -e $saclst or ! -e $sacch) {die(" No $saclst file\n");}
if (! -e $phtimes) {die("No $phtimes file\n");}
if (! -e $asc2sac) {die("No $asc2sac file\n");}


if ($opt_m) { # extract event information from CMTSOLUTION
  if (!-f $opt_m) {die(" No cmtsolution file $opt_m\n");}
#  ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,$hdur,$elat,$elon,$edep)=&get_cmt($opt_m);
#  ($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1)=&tdiff($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift);
($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1,$evid,$tshift,$hdur,
         $elat,$elon,$edep,undef,undef,undef,undef,undef,undef,undef,@evnm) = 
  get_cmt_long2($opt_m);
}

foreach $file (@ARGV) {

  print "\nProcessing file            $file\n"; 
  if (! -f $file) {die(" No such file : $file\n");}

  # transfer ascii file to sac file according to name of file
  if (not $opt_S) {# ascii file
    print "Transfer ascii file to sac file\n";
    system("$asc2sac $file >/dev/null");
    $file = $file.".sac";
  }
  (undef,$origin_time) = split(" ",`$saclst o f $file`);
  if (not defined $origin_time) {die("Check if the file is in SAC format");}

  ($filename) = split(" ",`basename $file`);
  ($sta,$net,$comp)=split(/\./,$filename);
  if (not $opt_d) {$outfile = $file.$ext;}
  else {$outfile = "$opt_d/${filename}${ext}";}
  if ($ext or $opt_d) {system("\\cp -f $file $outfile");}
  print "Output to file $outfile\n";

  if ($opt_m) {
    system("$sacch kevnm '@evnm' kuser0 $evid o 0 f $outfile >/dev/null"); }

 if ($opt_h) { # convolve source time function in SAC
    if (not defined $hdur) {die("Specify half duration for convolving source time function\n");} 
    print "Convolve source time function, hdur = $hdur\n";
    system("/opt/seismo-util/bin/convolve_stf g $hdur $outfile; mv -f $outfile.conv $outfile");}

  if ($opt_c) { open(SAC,"|sac2000");}else {open(SAC,"|sac2000 > /dev/null");}
  print SAC "echo on\n";
  print SAC "r $outfile\n";

  if ($opt_m) { # add event location and original time
    print "Add event information\n";
    if($opt_o){$tmp_o = -$opt_o; print SAC "ch allt $tmp_o\n";}
    print SAC "ch nzyear $oyear1 nzjday $ojday1 nzhour $ohr1 nzmin $omin1 nzsec $osec1 nzmsec $omsec1\n";
    print SAC "ch o 0\n";
    print SAC "ch evla $elat evlo $elon evdp $edep \n wh\n";}

  if ($opt_a) { # add station information
    print "Add station information\n";
    print SAC "ch kstnm $sta knetwk $net kcmpnm $comp\n ";
    if    ($comp=~/E/) {print SAC "ch cmpaz 90 cmpinc 90\n";}
    elsif ($comp=~/N/) {print SAC "ch cmpaz 0 cmpinc 90\n";}
    elsif ($comp=~/Z/) {print SAC "ch cmpaz 0 cmpinc 0 \n";}
    if (! -f $opt_a ) {die(" No such files: $opt_a");}
    ($sta_name,$sta_net,$sta_lat,$sta_lon,$sta_ele,$sta_bur)=split(" ",`egrep '$sta +$net' $opt_a`);
    if (not defined $sta_name) {
      print("No such station $sta+$net in $opt_a file\n");
      ($sta_name,$sta_net,$sta_lat,$sta_lon,$sta_ele,$sta_bur)=split(" ",`egrep '$sta' $opt_a`);
      if (not defined $sta_name) {die("No such station $sta_name in $opt_a file");}}
    print SAC "ch stla $sta_lat stlo $sta_lon stel $sta_ele stdp $sta_bur\nwh\n";}

  if ($opt_l) {
     print SAC "setbb begin ( max &1,b ( &1,o + $lmin ) ) \n";
     print SAC "setbb end   ( min &1,e ( &1,o + $lmax ) ) \n";
     print SAC "cut %begin% %end% \n";
     print SAC "r $outfile\n";
     print SAC "cut off\n";
 }

  print SAC "w over\nquit\n";
  close(SAC);

  if ($opt_c) { open(SAC,"|sac2000");}else {open(SAC,"|sac2000 > /dev/null");}
  print SAC "echo on\n r $outfile\n ";

  if ($opt_t ){ # filter record
    print "Filtering records...\n";
    print SAC "rtrend\n rmean\n taper\n";
    printf SAC ("bp n $opt_P p 2 co %10.5f %10.5f\n",$f1,$f2);
  }

  if ($opt_i) {# convolve with instrument response
    $pz= find_pz($opt_i,$outfile);
    print "Convolving instrument response using $pz...\n";
    if (not $pz) {next;}
    if (! -f $pz) {die("Not a pz file $pz\n");}
    print SAC "transfer from none to polezero s ${pz}\n";
    print SAC "w over\n";
  }

  if ($opt_p) { # pick arrival
    print SAC "w $outfile \nquit\n";
    close(SAC);
    if ($opt_c) { open(SAC,"|sac2000");}else {open(SAC,"|sac2000 > /dev/null");}
    print SAC "echo on\n r $outfile\n ";

    print "Picking arrivals...\n";
    (undef,$dist, $edep)=split(" ",`$saclst dist evdp f $outfile`);
    if (abs($dist-$undef)< $eps) {die("No dist is defined\n");}
    if (abs($edep-$undef) <$eps) {die("No event depth is defined, not able to calculate trav\n");}
    (undef,$tp) = split(" ",`trav.pl p /opt/seismo-util/source/fk1.0/socal $edep $dist`);
    if ($? != 0) {die("Error getting tp header for $outfile\n");}
    (undef,$ts) = split(" ",`trav.pl s /opt/seismo-util/source/fk1.0/socal $edep $dist`);
    if ($? != 0) {die("Error getting ts header for $outfile\n");}
    print SAC "evaluate to tmp1 $tp - $tshift\n";
    print SAC "evaluate to tmp2 $ts - $tshift\n";
    print SAC "ch t1 %tmp1% t2 %tmp2%\n";
    print SAC "ch kt1 P kt2 S\n wh\n";
  }

  if ($opt_y) { # add arrival time for surface waves
    @numbers=split(/\//,$opt_y); $npairs = floor(@numbers/2);
    if (@numbers != $npairs*2) {die("Enter -y in1/slope1/in2/slope2/...\n");}
    for ($i=0;$i<$npairs;$i++) {
      $int = $numbers[$i*2]; $slope = $numbers[$i*2+1];
      print "Add arrival time for waves with group velocity: $slope\n";
      (undef,$dist,$begin,$end)=split(" ",`$saclst dist b e f $outfile`);
      if (abs($dist-$undef) < $eps) {die("Not defined dist\n");}
      $h1 = 3+$i;  $h1 = "t$h1"; $k1 = "k$h1"; $v1 = $int + $dist/$slope;
      (undef,$v1,undef) = sort {$a <=> $b} ($begin, $v1 ,$end);
      printf SAC ("ch $h1 %12.2f\n",$v1);
      print SAC "ch $k1 $h1\n";
    }
  }
  if ($dt) {print SAC "interp delta $dt\n";}
  print SAC "w over\nquit\n";
  close(SAC);

  if ($opt_w) {
    if (not $opt_m) {die("Specify CMTSOLUTION for -w option\n");}
    ($mrr,$mtt,$mpp,$mrt,$mrp,$mtp,$exp) = get_moment($opt_m);
    system("$sacch resp0 $mrr resp1 $mtt resp2 $mpp resp3 $mrt resp4 $mrp resp5 $mtp isynth $exp f $outfile >& /dev/null");
  }

}
print " \n     DONE\n\n";
     
#**********************************************************

sub find_pz {
  my($pz_dir,$oldfile)=@_;
  my ($i);
  my($sta,$net,$comp)=split(/\./,`basename $oldfile`);
  my($saclst) = "/opt/seismo-util/bin/saclst";
  my(undef,$year,$jday,$hour,$min,$sec)=split(" ",`$saclst nzyear f $oldfile`);
  my($mon,$day)=&jday2mday($year,$jday);
  if ($mon<10) {$mon="0".$mon;}
  if ($day<10) {$day="0".$day;}
  my($pzdate)=$year.$mon.$day;
  my(@pzfiles)=`ls -1 $pz_dir | egrep -i "${net}_${sta}_${comp}_*" `;
  if (@pzfiles ==0) {
    print ("No pz file for station $sta, skipped\n");
    return("");}
  for ($i=0;$i<@pzfiles;$i++) {
    chomp($pzfiles[$i]);
    my(undef,undef,undef,$start,$end)=split("_",$pzfiles[$i]);
    if ($pzdate ge $start and $pzdate le $end) {$pzfile = $pzfiles[$i];last;}
    if ($i == @pzfiles-1) {die("No proper pzfiles\n");}
  }
  $pzfile=$pz_dir."/".$pzfile;
  return($pzfile);
}

sub jday2mday {
  my($year,$jday)=@_;
  $time = timegm(1,0,0,1,0,$year);
  $time=$time+($jday-1)*24*60*60;
  my @gmt_time=gmtime($time);
  return($gmt_time[4]+1,$gmt_time[3]);
}

sub get_moment {
  my($cmt_file) = @_;
  open(CMT_SUB,"$cmt_file") || die("Can not open $cmt_file");
  my(@cmt) = <CMT_SUB>; close(CMT_SUB);
  my(@m) = (); my(@exp) = ();
  my($i);
  for ($i=7; $i<13; $i++) {
    my(undef, $value) = split(" ",$cmt[$i]);
    my($val,$expon) = split(/\+/,$value);
    chop($val);
    $exp[$i-7] = $expon;
    $m[$i-7] = $val;
  }
  ($max) = sort {$b <=> $a} @exp;
  for ($i=0;$i<6; $i++) {
    my($exp1) = $exp[$i] - $max;
    $m[$i] = $m[$i] * "1.0e$exp1";
  }
  return(@m[0..5],$max);
}
