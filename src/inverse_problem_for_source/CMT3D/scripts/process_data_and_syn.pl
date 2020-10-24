#!/usr/bin/perl -w

# again this is partly based on Carl's main script
use Getopt::Std;
use POSIX;
use vars qw($opt_0 $opt_1 $opt_2 $opt_3 $opt_t $opt_x);
use List::Util qw[min max];

sub Usage{
   print STDERR <<END;

     process_data_and_syn.pl
        -x -d data_dir,data_ext -s syn_dir,syn_ext
        -0 -m cmt_file -a sta_file -n sps
        -1 -b bmin(-21)
        -2 -t Tmin/Tmax
        -3 [-t Tmin/Tmax]


    Note: -h is used to convolve the raw synthetics with triangular stf
          (taken out for now)
          -0: pre-processing (cmt,sta,pick)
          -1: cutting and padding (bmin)
          -2: filtering
          -3: flexwin_input
END
   exit(1);
 }

if (@ARGV == 0) {Usage();}
# cmt_file, sta_file, data_dir, data_ext, syn_dir, syn_ext, Tmin, Tmax, opt_i, nsps, bmin (-20)
if (!getopts('m:a:d:s:01b:23t:n:x')) {die("Check input options\n");}

if ($opt_0) {$pre_process = 1;} else {$pre_process = 0;}
if ($opt_1) {$cut = 1;} else {$cut = 0;}
if ($opt_2) {$band_pass = 1;} else {$band_pass = 0;}
if ($opt_b) {$bmin = $opt_b;} else {$bmin = 200;} # some random number
if ($opt_3) {$flexwin_input=1;} else {$flexwin_input = 0;}

# data dir and extension
if ($opt_d) {
  ($data_dir,$data_ext) = split(/\,/,$opt_d);
  if ($data_ext) {$data_ext = ".$data_ext";}
  else {$data_ext="";}}
if ($opt_s) {
  ($syn_dir,$syn_ext) = split(/\,/,$opt_s);
  if ($syn_ext) {$syn_ext=".$syn_ext";}
  else {$syn_ext="";}}

# check CMTs and STATIONS file
if ($pre_process) {
  if ($opt_m) {$cmt_file=$opt_m;} else {die("-m option\n");}
  if ($opt_a) {$sta_file=$opt_a;} else {die("-a option\n");}
  if (not -f $cmt_file or not -f $sta_file) {die("Check $cmt_file or $sta_file\n");}
  if ($opt_n) {$nsps="-s $opt_n";} else {$nsps = "-s 20";}

# uncomment this if you choose to convolve hdur here
#  (undef,undef, $hdur) = split(" ",`grep half $cmt_file`);
#  if (abs($hdur)<0.1) {die("Check half duration for -h option in syn processing\n");}
}

if ($opt_d) {$dirdat_pp= "${data_dir}_PP";}
if ($opt_s) {$dirsyn_pp = "${syn_dir}_PP";}

# period range
if ($band_pass or $flexwin_input) {
  if (not $opt_t) {die("Input -tTmin/Tmax for bandpassing\n");}
  ($Tmin,$Tmax) = split("/",$opt_t);
  $Trange="$Tmin/$Tmax";
  $sTmin = sprintf("T%3.3i",$Tmin);
  $sTmax = sprintf("T%3.3i",$Tmax);
  $Ttag = "${sTmin}_${sTmax}";
  if ($opt_d) {$dirdat_fil = "${data_dir}_${Ttag}";}
  if ($opt_s) {$dirsyn_fil = "${syn_dir}_${Ttag}";}}

if ($flexwin_input) {
  $flex_in="flexwin_${Ttag}.input";
  $flex_out="flexwin_${Ttag}.stdout";
  $flex_dir="flexwin_${Ttag}";}

if ($cut or $flexwin_input) {
    if (not $opt_s or not $opt_d) {die("Cutting (-1): data_dir and syn_dir have to be both present\n");}}

# write the C-shell script to file
$cshfile = "process_ds.csh";
print "\nWriting to $cshfile ...\n";
open(CSH,">$cshfile");

# ---------------------------------------------------------------
if ($pre_process) { # 0: pre-processing (data_dir or syn_dir)
  print CSH "echo 'Pre-processing data and/or syn ...'\n";
  if ($opt_s and -d $syn_dir) {
    if (-d $dirsyn_pp) {print "**** $dirsyn_pp already exists *****, replacing ...\n";}
    print CSH "process_trinet_syn_new.pl -S -m $cmt_file -a $sta_file $nsps -p -d $dirsyn_pp ${syn_dir}/*H?*${syn_ext} \n";}

  if ($opt_d and -d $data_dir) {
    if (-d $dirdat_pp) {print "**** $dirdat_pp already exists *****, replacing ...\n";}
    print CSH "process_cal_data.pl -m $cmt_file $nsps -p -d $dirdat_pp ${data_dir}/*H?${data_ext} \n";}
}

# ----------------------------------------------------------------
if ($cut) { # 1: padding and cutting
  # pad zeros
  if ($opt_b) {
    print CSH "echo 'Padding zeros to data and syn ...'\n";
    print CSH "pad_zeros.pl -s -l $bmin ${dirdat_pp}/*H?$data_ext >/dev/null\n";
    print CSH "pad_zeros.pl -s -l $bmin ${dirsyn_pp}/*H?*$syn_ext >/dev/null\n";}

  $cutdat_file="$dirdat_pp/CUT_DATA_PP";
  $cutsyn_file="$dirsyn_pp/CUT_SYN_PP";
  open(CUTDAT,">$cutdat_file");
  open(CUTSYN,">$cutsyn_file");

  # matching data and syn
  @files = glob("${dirdat_pp}/*H?$data_ext");
  $nfile = @files;
  print CSH "echo \"Cutting ... $nfile data files to line up with synthetics\"\n";
  print CSH "sac <<EOF\n";

  foreach $datfile (@files) {

    (undef,$net,$sta,$chan) = split(" ",`saclst knetwk kstnm kcmpnm f $datfile`);
    ($comp) = split(" ",`echo $chan | awk '{print substr(\$1,3,1)}'`);

    $synfile = "${dirsyn_pp}/${sta}.${net}.BH${comp}${syn_ext}";

    if (not -f $synfile) {
      print "**** Missing synthetic pair $synfile for $datfile ***, skipping ...\n";
    } else {
#      print "$datfile $synfile\n";

      # get info on data and synthetics
      (undef,$bd,$ed,$deltad) = split(" ",`saclst b e delta f $datfile`);
      (undef,$bs,$es,$deltas) = split(" ",`saclst b e delta f $synfile`);

      # dt should be the same for both records ALREADY
      if (abs($deltad-$deltas) > 0.00001) {
  print "DT values are not close enough: $deltad, $deltas\n";
  die("fix the DT values\n");
      } else {$dt =$deltad;}
      $b = min(max($bd,$bs),$bmin+1);
      $e0 = min($ed,$es);
      $tlen = $e0 - $b;
#      $tlen = $tfac * $tlen0;
      $e = $b + $tlen;
      $npt = int( ($e-$b)/$dt );

      @synfiles = glob("${dirsyn_pp}/${sta}.${net}.BH${comp}*${syn_ext}");
      print CUTDAT "$datfile $b $e $npt $dt\n";
      print CUTSYN "$synfile $b $e $npt $dt\n";

      print CSH "cut $b n $npt\n";
      print CSH "r $datfile @synfiles\n";
      print CSH "cut off\n";
      print CSH "w over \n";
    }
  }
  close(CUTDAT);
  close(CUTSYN);
  print CSH "quit\n";
  print CSH "EOF\n";
  print CSH "rm -f tmpPadZeroToSyn.sac\n";
}
#----------------------------------------------------------------------
if ($band_pass) {
  print CSH "echo 'Bandpassing data and/or synthetics ...'\n";
  if ($opt_s) {
    print CSH "process_trinet_syn_new.pl -S -t $Trange -d $dirsyn_fil $dirsyn_pp/*.?H?*$syn_ext \n";}
    print CSH "rotate.pl $dirsyn_fil/*.?H[E1]*${syn_ext}\n";

  if ($opt_d) {
    print CSH "process_cal_data.pl -i none -t $Trange -d $dirdat_fil -x d $dirdat_pp/*.?H?$data_ext \n";}
    print CSH "rotate.pl $dirdat_fil/*.?HE${data_ext}.d\n";
}
# remaining: record the failed ones from output

#----------------------------------------------------------------------
close(CSH);
print "\n";
system("chmod a+x $cshfile");
if ($opt_x) {system("$cshfile");}

if ($flexwin_input) {
  if (not $opt_x) {die("Make sure c-shell script run before flexwin_input\n");}

  open (FLX,">temp.input");
  # matching data and syn (ZRT)
  @files = glob("${dirdat_fil}/*H[ZRT]${data_ext}.d");
  $nfile = @files;  $nn = 0;
  foreach $datfile (@files) {

    (undef,$net,$sta,$chan) = split(" ",`saclst knetwk kstnm kcmpnm f $datfile`);
    ($comp) = split(" ",`echo $chan | awk '{print substr(\$1,3,1)}'`);
    $synfile = "${dirsyn_fil}/${sta}.${net}.BH${comp}${syn_ext}";
    if (not -f $synfile) {
      print "**** Missing synthetic pair $synfile for $datfile ***, skipping ...\n";
    } else {
      print FLX "$datfile\n$synfile\n$flex_dir/${sta}.${net}.BH${comp}.${Ttag}\n";
      $nn++;
    }
  }
  system("echo $nn > $flex_in; cat temp.input >> $flex_in; rm -f temp.input");
  system("mkdir -p $flex_dir");
  system("cp flexwin_user_fun/PAR_FILE_$Ttag PAR_FILE");
  system("cp flexwin_user_fun/user_functions.f90 ../flexwin; cd ../flexwin; make");
  if (not -f "iasp91.hed") {system("ln -s ../flexwin/ttimes_mod/iasp91.hed");}
  if (not -f "iasp91.tbl") {system("ln -s ../flexwin/ttimes_mod/iasp91.tbl");}
  print "Running flexwin with input file $flex_in ...\n";
  system("flexwin < $flex_in > ${flex_out}");
}

