#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;

sub Usage{
print STDERR <<END;

Usage:   process_syn.pl
     -S -m CMTFILE -h -o offset -lStart/End -tTmin/Tmax -f -P n/p
     -i Dir -A amp_factor -p -a STAFILE -s sps -y t0/V0/t1/V1...
     -c -d OutDir -x Ext  synthetics_files
where
    -S already sac file, will not apply ascii2sac.csh
    -m CMTSOLUTION file for event origin, location and half duration
    -h convolve a triangle stf with half duration from -m
    -o moves the synthetics back by 'offset' from centroid time.
    -l lmin/lmax   start/end cut of trace
    -t Tmin/Tmax -- shortest/longest period of bandpass filter
    -f -- apply transfer filter instead of butterworth filter for -t
    -P -- number of poles and passes in butterworth filter(default 4/2)
    -i Dir -- convolve synthetics with instrument response in Dir
    -A amp_factor -- multiple all files with an amplification factor
    -p -- pick P and S arrivals according to iasp91 model to t1/t2 headers
    -a STAFILE -- add station information from STAFILE, use "-a none"
       to get information from Vala''s station file
    -s sps -- resampling rate for synthetics (default is 1.0)
    -y t0/V0/t1/V1 -- add to t3/t4 the arrival times of t + dist/V
    -c -- show sac output
    -d OutDir -- directory to dump output (default: same as input files)
    -x Ext -- add extension Ext

    names of files -- name of syn files to be processed

 NOTE: Please make sure that SAC, saclst and IASP91 packages are installed properly on 
       your system, and that all the environment variables are set properly before
       running the script.

    
    Qinya Liu, originally written in Oct 2002, most recently updated in Feb 2010.

END
exit(1);
}

@ARGV > 1 or Usage();

if (!getopts('Sm:ho:l:t:fi:pP:a:cd:x:vs:y:A:')) {die(" check input arguments\n");}

if ($opt_t) {($tmin, $tmax) = split(/\//,$opt_t);
             $f1 = 1./$tmax;  $f2=1./$tmin;}
if ($opt_x) {$ext=".".$opt_x;} else {$ext="";}
if (!$opt_s) {$dt=1.0;} else {$dt = 1.0/$opt_s;}
if (!$opt_P) {$poles=4;$pass=2;}
else{($poles,$pass)=split(/\//,$opt_P);
     if(not defined $pass or $pass<1){$pass=2;}}
if ($opt_l) {($lmin,$lmax) = split(/\//,$opt_l);} 
else {$lmin = 0; $lmax = 3600;}
if ($opt_a and not -f $opt_a) {$opt_a="/opt/seismo/data/STATIONS";}
if ($opt_o and not $opt_m) {die("Specify centroid time first\n");}
if ($opt_d and not -d $opt_d) {die("No such directory as $opt_d\n");}
if ($opt_i and not -d $opt_i) {die("No such directory as $opt_i\n");}
if ($opt_A) {if ($opt_A !~ /^\d/) {die("-A option should be numbers\n");}}
if ($opt_f) {
  if (not $opt_t) {die("-t tmin/tmax has to be specified for -f\n");}
  else {$f0=$f1*0.8;  $f3=$f2*1.2;}}

$undef = -12345;
$eps = 0.1;

$saclst="saclst";
$phtimes="./phtimes.csh";
$asc2sac="ascii2sac.csh";

#if (! -e $saclst)  {die(" No $saclst file\n");}
#if (! -e $phtimes) {die("No $phtimes file\n");}
#if (! -e $asc2sac) {die("No $asc2sac file\n");}

$min_hdur=1.0;

if ($opt_m) { # extract event information from CMTSOLUTION
  if (!-f $opt_m) {die(" No cmtsolution file $opt_m\n");}
  ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,$hdur,$elat,$elon,$edep) = get_cmt($opt_m);
  ($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1) = tdiff($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift);}

$sta_text="";

foreach $file (@ARGV) {

  print "\nProcessing file            $file\n"; 
  if (! -f $file) {die(" No such file : $file\n");}

  # transfer ascii file to sac file according to name of files
  if (not $opt_S) {# ascii file
    print "Transfer ascii file to sac file\n";
    system("$asc2sac $file >/dev/null");
    $file = $file.".sac";}

  (undef,$begin_time) = split(" ",`$saclst b f $file`);
  if (not defined $begin_time) {die("Check if the file is SAC format\n")};

  ($filename) = split(" ",`basename $file`);
  ($sta,$net,$comp)=split(/\./,$filename);
  if (not $opt_d) {$outfile = $file.$ext;}
  else {$outfile = "$opt_d/${filename}${ext}";}
  if ($ext or $opt_d) {system("\\cp -f $file $outfile");}
  print "Output to file $outfile\n";

  if ($opt_c) { open(SAC,"|sac");}else {open(SAC,"|sac > /dev/null");}
  print SAC "echo on\n";

  print SAC "r $outfile\n";
  if ($opt_m) { # add event location and original time
    print "Add event information\n";
    if($opt_o){$tmp_o = -$opt_o; print SAC "ch allt $tmp_o\n";}
    print SAC "ch nzyear $oyear1 nzjday $ojday1 nzhour $ohr1 nzmin $omin1 nzsec $osec1 nzmsec $omsec1\n";
    print SAC "ch o 0\n";
    print SAC "ch evla $elat evlo $elon evdp $edep \n wh\n";
    print SAC "w over \nquit\n";
    close(SAC);
    if ($opt_c) {open(SAC,"|sac");}else {open(SAC,"|sac > /dev/null");}
    print SAC "echo on\n r $outfile\n";}

  if ($opt_h and $hdur > $min_hdur) { # convolve source time function
    print SAC "quit\n";
    close(SAC);
    system("/opt/seismo-utils/bin/convolve_stf g $hdur $outfile");
    system("mv ${outfile}.conv ${outfile}");
    if ($opt_c) {open(SAC,"|sac");}else {open(SAC,"|sac > /dev/null");}
    print SAC "echo on\n r $outfile\n";}


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
      if (not defined $sta_name) {
         print " No such station as $sta in the $opt_a file\n";
         $sta_text .="$sta ";}}
    print SAC "ch stla $sta_lat stlo $sta_lon stel $sta_ele stdp $sta_bur\nwh\n";}

   print SAC "w over\nquit\n";
   close(SAC);
   if ($opt_c) { open(SAC,"|sac");}else {open(SAC,"|sac > /dev/null");}
   print SAC "echo on\n r $outfile\n ";

  if ($opt_l) { #cut record
     print SAC "setbb begin ( max &1,b ( &1,o + $lmin ) ) \n";
     print SAC "setbb end   ( min &1,e ( &1,o + $lmax ) ) \n";
     print SAC "cut %begin% %end% \n";
     print SAC "r $outfile\n";
     print SAC "cut off\n";}

  if ($opt_A) {print SAC "mul $opt_A\n";}

  if ($opt_t)  { # filter record
    print "Filtering records...\n";
    print SAC "rtrend\n rmean\n taper\n";
    if ($opt_f) {
      printf SAC ("trans from none to none freq %12.6f %12.6f %12.6f %12.6f\n",$f0,$f1,$f2,$f3);}
    else {
      printf SAC ("bp n $poles p $pass co %10.5f %10.5f\n",$f1,$f2);}
    print SAC "rtrend\n rmean\n taper\n";}
 
  if ($opt_i) {# convolve with instrument response
    print "Convolving instrument response...\n";
    $pzfile="SAC_PZs_${net}_${sta}_${comp}_";
    @nfiles=`ls -l $opt_i/${pzfile}* | awk '{print \$9}'`;
    if (@nfiles != 1) {die("Pzfile is incorrect: \n@nfiles files\n");}
    $pz=$nfiles[0]; chomp($pz);
    if (! -f $pz) {die("Not a pz file $pz\n");}
    print SAC "transfer from none to polezero s ${pz}\n";
    print SAC "w over\n";}

  if ($opt_s) {print SAC "interp delta $dt\n";}

  if ($opt_p) { # pick arrival
    print "Picking arrivals...\n";
    (undef,$gcarc)=split(" ",`$saclst gcarc f $outfile`);
    if (abs($gcarc-$undef)< $eps) {die("No gcarc is defined\n");}
    ($Pph,$Ptime)=split(" ",`$phtimes $edep $gcarc P`);
    ($Sph,$Stime)=split(" ",`$phtimes $edep $gcarc S`);
    print SAC "evaluate to tmp1 $Ptime - $tshift\n";
    print SAC "evaluate to tmp2 $Stime - $tshift\n";
    print SAC "ch t1 %tmp1% t2 %tmp2%\n";
    print SAC "ch kt1 $Pph kt2 $Sph\n wh\n";}

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
      print SAC "ch $k1 $h1\n";}}

  print SAC "w over\nquit\n";
  close(SAC);

}
if ($sta_text) {
   print "\n\nNotice the following stations are not found in the station file\n   $sta_text\n\n";}
print "DONE\n";


#****************************************************************
sub mday2jday {
    my($oyear,$omonth,$oday)=@_;
    $omonth = $omonth-1;  #months range from 0..11
    $time_sec = timegm(3,3,3,$oday,$omonth,$oyear);
    @t = gmtime($time_sec);
    $jday = @t[7];
    $jday += 1;
    return ($jday);
}

sub tdiff {
    my ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tadd)=@_;
    #  die("$oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tadd\n");
    $time = timegm($osec, $omin, $ohr, $oday , $omonth-1, $oyear);
    $time += ($tadd +$omsec/1000); #event_time in machine format
    $msec = sprintf("%03.0f",($time - (floor($time)))*1000);
    $time = floor($time);
    #event-time:
    ($sec, $min, $hr, $day , $month, $year,$weekday,$jday) = gmtime($time);
    $month += 1;
    $year += 1900;
    $jday +=1;
    return ($year,$jday,$hr,$min,$sec,$msec);
}

sub get_cmt {
    local ($cmt_file)=@_;
    open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
    @cmt = <CMT>;
    local($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1)=split(" ",$cmt[0]);
    local($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
    local(undef,undef,$tshift)=split(" ",$cmt[2]);
    local(undef,undef,$hdur)=split(" ",$cmt[3]);
    local(undef,$elat)=split(" ",$cmt[4]);
    local(undef,$elon)=split(" ",$cmt[5]);
    local(undef,$edep)=split(" ",$cmt[6]);
    local(undef,$Mrr)=split(" ",$cmt[7]);
    local(undef,$Mtt)=split(" ",$cmt[8]);
    local(undef,$Mpp)=split(" ",$cmt[9]);
    local(undef,$Mrt)=split(" ",$cmt[10]);
    local(undef,$Mrp)=split(" ",$cmt[11]);
    local(undef,$Mtp)=split(" ",$cmt[12]);
    close(CMT);
    return ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift,$hdur,$elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp);
}

