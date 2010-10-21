#!/usr/bin/perl

use Time::Local;
use Math::Trig;
use Getopt::Std;
use POSIX;
use lib '/opt/seismo-util/lib/perl/';
use CMT_TOOLS;

sub Usage{

  print STDERR <<END;

  process_trinet_data_new.pl  -m CMTFILE -l Start/End -t Ts/Tl 
           -[i|I] Dir -p -x Ext -d OutDir -s Sps -P n
           -c -y t1/s1/...  -v num_of_zeros names_of_files
  with

       -m -- Using CMTSOLUTION file to set event info
       -l Start/End of trace to be cut
       -t Ts/Tl -- shortest and longest period
       -[i|I] Dir-- transfer recordings to displacement with pz files in Dir
        (use none to set to defaulted Jascha''s dir)
        (use -I instead to only perform transfer no bp)
       -p -- pick event p and s arrival in sac (use socal model)
       -x Ext -- extension of new file name
       -d Dir -- directory to put the output files
       -s Sps -- sampling rate per second. (default 20)
       -P n -- number of poles for the butterworth filter (default 4)
       -c  -- check sac output on the screen
       -y t0/slope/t0/slope... (t3/t4/...)(for SC, try -20/0.15/11-/0.2)
       -v # -- convert the data to either velocity or displacement by specifying
               the number of zeros [default no change from pz files]

       name_of_files --- name of sac files to be processed

Notice: 1) main difference from process_data_new.pl is to specify 
        the pz files directory and the default sampling rate is for BH component.
        2) We always use the polezeros files corresponding to the correct
        station, component in the pz directory. Polezero files are in the
        format of StationName_CompName_StartTime_EndTime and recommended
        pz directory is /opt/seismo-util/data/TriNet_PZs/, in other
        cases, it would not work. Give "none" for this default directory.
        Output data files are in the units of meters.
        4) please ignore any warning saying Number 4003 too large, it is
        an intrinsic bug of the sac version we use.

END
exit(1)
}

if (@ARGV == 0) { Usage(); }

if (!getopts('m:l:t:i:I:px:d:s:P:cy:v:')) {die('Check input arguments\n');}
if ($opt_t) {($tmin, $tmax) = split(/\//,$opt_t);$f1 = 1./$tmax;  $f2=1./$tmin;}
if ($opt_d) {$out_dir=$opt_d;if (not -d $out_dir) {mkdir($out_dir,0777);}}
if ($opt_x) {$ext='.'.$opt_x;}else {$ext="";}
if (!$opt_s) {$dt=1.0/20;} else {$dt = sprintf("%2.2e",1.0/$opt_s);}
if (!$opt_P) {$opt_P=4;}
if ($opt_l) {($lmin,$lmax) = split(/\//,$opt_l);} else {$lmin = 0; $lmax = 200;}
if ($opt_i and ! -d $opt_i) {$opt_i="/opt/seismo-util/data/TriNet_PZs";}
if ($opt_i and $opt_I) {die(" -I and -i can't coexist\n");}
if ($opt_I and ! -d $opt_I) {$opt_i="/opt/seismo-util/data/TriNet_PZs";}
if (not $opt_v) {$opt_v = 0;}

$saclst="/opt/seismo-util/source/saclst/saclst";
$phtimes="/opt/seismo-util/bin/phtimes";
if (! -e $saclst) {die(" No $saclst file\n");}
if (! -e $phtimes) {die("No $phtimes file\n");}

$eps=1e-5; $undef=-12345.0;

foreach $file (@ARGV){
  if (! -e $file) {die("No such file to be processed!!\n");}
  print "Processing data file $file \n";
  if ( ! $opt_d ) {$outfile = $file.$ext;}
  else { ($filename) = split(" ",`basename $file`);
	 $outfile = "$out_dir/${filename}${ext}";}
  if ($ext or $opt_d) {system("\\cp -f $file $outfile");}

 # process data
  if ($opt_c) {open(SAC,"|sac2000  ") || die("Can't open sac2000\n");}
  else {open(SAC,"|sac2000 >sac_out.log ") || die("Can't open sac2000\n");}
  print SAC "echo on\n";
  print SAC "r $outfile\n";
  
  if ($opt_m) { # USE CMTSOLUTION FILE
    #   ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,undef,$elat,$elon,$edep)=&get_cmt($opt_m);
    #   ($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1)=&tdiff($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift);
    ($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1,$evid,$tshift,$hdur,
     $elat,$elon,$edep,undef,undef,undef,undef,undef,undef,undef,@evnm) =
       get_cmt_long2($opt_m);
    print SAC "ch o gmt $oyear1 $ojday1 $ohr1 $omin1 $osec1 $omsec1\n";
    print SAC "evaluate to  tmp 0 - &1,o\n";
    print SAC "ch allt %tmp% iztype io\n";
    print SAC "ch evla $elat evlo $elon evdp $edep \n" ;
    print SAC "w $outfile \nquit\n";
    close(SAC);
    system("/opt/seismo-util/bin/sacch kevnm '@evnm Data' kuser0 $evid f $outfile > /dev/null");
   if ($opt_c) {open(SAC,"|sac2000 ") || die("Can't open sac2000\n");}
   else {open(SAC,"|sac2000 > sac_out.log") || die("Can't open sac2000\n");}
   print SAC "echo on\n";
   print SAC "r $outfile\n";

 }

  if ($opt_l){  # cut record and rtrend and rmean
    print "    Cut record from o+$lmin to o+$lmax\n";
    (undef,$tmp_o)=split(" ",`$saclst o f $outfile`);
    if (abs($tmp_o-$undef) < $eps) {die("Undefined o, can't cut record\n");}
     print SAC "setbb begin ( max &1,b ( &1,o + $lmin ) ) \n";
     print SAC "setbb end   ( min &1,e ( &1,o + $lmax ) ) \n";
     print SAC "cut %begin% %end% \n";
     print SAC "r $outfile\n cut off \n w over \n";
     print SAC "rmean\n rtrend\n taper\n";
  }

  if ($opt_i or $opt_I) {  # transfer instrument response
    if (! $opt_t) {die(" Specify bandpass range by -t\n");}
    $f0=$f1*0.8;
    $f3=$f2*1.2;
    if ($opt_i) {$pz_dir = $opt_i;} else {$pz_dir = $opt_I;}
    $pzfile= find_pz($pz_dir,$file,$opt_v);
    if (not $pzfile) {
      if (!$opt_d and !$ext) {die("Check $pz_dir for $file\n");}
      else {
	print "Check $pz_dir for polezero file of $file-- skip\n";
	system("\\rm -f $outfile");  next;}}
    print "    Transfer instrument response using pz file $pzfile\n";
    print SAC "rmean \n rtrend \n taper\n";
    printf SAC ("trans from polezero s %25s to none freq %12.6f%12.6f%12.6f%12.6f \n",$pzfile,$f0,$f1,$f2,$f3);
  }

  if ($opt_t and not $opt_I) {# filter records
    print "    Filter record at $tmin and $tmax second\n";
    printf SAC ("bp n %4d p 2 cor %12.6f %12.6f\n",$opt_P,$f1,$f2);
}

  
  if ($opt_p) {   # add p and s arrival info
    print SAC "w $outfile \nquit\n";
    close(SAC);
    if ($opt_c) {open(SAC,"|sac2000  ") || die("Can't open sac2000\n");}
    else {open(SAC,"|sac2000 >sac_out.log ") || die("Can't open sac2000\n");}
    print SAC "echo on\n";
    print SAC "r $outfile\n";
    print "Picking arrivals...\n";
    (undef,$dist,$edep)=split(" ",`$saclst dist evdp f $outfile`);
    if ($dist eq "" or $edep eq "") {die("Error taking dist and depth header \n");}
    if (abs($dist-$undef)< $eps) {die("No dist is defined\n");}
    if (abs($edep-$undef)< $eps) {die("No event depth is defined, not able to calculate trav\n");}
    (undef,$tp) = split(" ",`trav.pl p /opt/seismo-util/source/fk1.0/socal $edep $dist`);
    if ($? != 0) {die("Error getting tp header for $outfile\n");}
    (undef,$ts) = split(" ",`trav.pl s /opt/seismo-util/source/fk1.0/socal $edep $dist`);
    if ($? != 0) {die("Error getting tp header for $outfile\n");}
    if (not $opt_m) {$tshift = 0;}
    print SAC "evaluate to tmp1 $tp - $tshift\n";
    print SAC "evaluate to tmp2 $ts - $tshift\n";
    print SAC "ch t1 %tmp1% t2 %tmp2%\n";
    print SAC "ch kt1 P kt2 S\n";
    print SAC "wh\n";
  }

   # add arrival time
  if ($opt_y) {
    @numbers=split(/\//,$opt_y); $npairs = floor(@numbers/2);
    if (@numbers != $npairs*2) {die("Enter -y in1/slope1/in2/slope2/...\n");}
    for ($i=0;$i<$npairs;$i++) {
      $int = $numbers[$i*2]; $slope = $numbers[$i*2+1];
      print "   Add arrival time for waves with group velocity: $slope\n";
      (undef,$dist,$begin,$end)=split(" ",`$saclst dist b e f $outfile`);
      if (abs($dist-$undef) < $eps) {die("Not defined dist\n");}
      $h1 = 3+$i;  $h1 = "t$h1"; $k1 = "k$h1"; $v1 = $int + $dist/$slope;
      (undef,$v1,undef) = sort {$a <=> $b} ($begin, $v1 ,$end);
      printf SAC ("ch $h1 %12.2f\n",$v1);
      print SAC "ch $k1 $h1\n";
    }
  }

if ($opt_s){  print SAC "interp delta $dt\n";}
  print SAC "w $outfile\n";
  print "   write file $outfile\n";
  print SAC "echo off\nquit\n";
  close(SAC);
  system("\\rm -f sac_out.log");
}


#*************************************************************


sub find_pz {
  my($pz_dir,$oldfile,$opt_v)=@_;
  my ($i);
  my(undef,$sta,$net,$comp)=split(" ",`$saclst kstnm knetwk kcmpnm f $oldfile`);
  my($saclst) = "/opt/seismo-util/bin/saclst";
  my(undef,$year,$jday,$hour,$min,$sec)=split(" ",`$saclst nzyear f $oldfile`);
  my($mon,$day)=&jday2mday($year,$jday);
  if ($mon<10) {$mon="0".$mon;}
  if ($day<10) {$day="0".$day;}
  my($pzdate)=$year.$mon.$day;
  my(@pzfiles)=`ls -1 $pz_dir | egrep -i "${net}_${sta}_${comp}_*" `;
  if (@pzfiles ==0) { 
    print ("No pz file for station $sta: ls -1 $pz_dir | egrep -i \"${net}_${sta}_${comp}_*\" \n"); 
    return("");}
  for ($i=0;$i<@pzfiles;$i++) {
    chomp($pzfiles[$i]);
    my(undef,undef,undef,$start,$end)=split("_",$pzfiles[$i]);
    if ($pzdate ge $start and $pzdate le $end) {$pzfile = $pzfiles[$i];last;}
    if ($i == @pzfiles-1) {die("No proper pzfiles\n");}
  }
  $pzfile=$pz_dir."/".$pzfile;
  if ($opt_v and $opt_v != 2 and $opt_v != 4) {die("Check if -v has the right value [2/4]\n");}
  if ($opt_v) {system("sed '/ZEROS/s/3/$opt_v/' $pzfile > tmp_pzfile");
  $pzfile = "tmp_pzfile";}
  return($pzfile);
}

sub jday2mday {
  my($year,$jday)=@_;
  $time = timegm(1,0,0,1,0,$year);
  $time=$time+($jday-1)*24*60*60;
  my @gmt_time=gmtime($time);
  return($gmt_time[4]+1,$gmt_time[3]);
}
