package CMT_TOOLS;

#use warnings;
use Exporter;
use Time::Local;
use POSIX;

@ISA = ('Exporter');
@EXPORT = ('get_cmt','get_cmt_long','get_cmt_short','get_cmt_moment',
	   'get_cmt_long2', 'get_cmt_evnm', 'get_cmt_hdur_evdp',
	   'get_cmt_location','cmt2dc', 'cmt2dc_simple','moment2dc',
           'tdiff', 'tdiff2', 'mday2jday', 'jday2mday'
	  );

$BIN = "/opt/seismo-util/bin";

# get_cmt : all your want from a cmt_file except event name
sub get_cmt {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
  my(undef,undef,$evid) = split(" ",$cmt[1]);
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);

  return($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$evid,$tshift,$hdur,
         $elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp);

}


# get_cmt_long  : all you want from cmt_files including event name
sub get_cmt_long {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1,undef,undef,undef,$mw,undef,@evnm)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
  my(undef,undef,$evid) = split(" ",$cmt[1]);
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);

  return($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$evid,$tshift,$hdur,
         $elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp,$mw,@evnm);

}

# get_cmt_long2 : same as get_cmt_long except that the time shift has been
#                 added to the output reference time.
sub get_cmt_long2 {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1,undef,undef,undef,$mw,undef,@evnm)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec="0.$omsec";
  $omsec = floor($omsec * 1000);
  my(undef,undef,$evid) = split(" ",$cmt[1]);
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);

  ($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1)=&tdiff($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift);

  return($oyear1,$ojday1,$ohr1,$omin1,$osec1,$omsec1,$evid,$tshift,$hdur,
         $elat,$elon,$edep,$Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp,$mw,@evnm);

}

# get_cmt_evnm  : all you want from cmt_files including event name
sub get_cmt_evnm {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my($pde,$oyear,$omonth,$oday,$ohr,$omin,$osec1,undef,undef,undef,$mw,undef,@evnm)=split(" ",$cmt[0]);
  my($osec,$omsec)=split(/\./,$osec1); $omsec=$omsec*10;
  my(undef,undef,$evid) = split(" ",$cmt[1]);
  my(undef,undef,$tshift)=split(" ",$cmt[2]);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);
  my($evnm) = join(" ",@evnm);
  return($evnm);

}


# get_cmt_short : simply outputs M and location
sub get_cmt_short {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);

  return($Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp,$elat,$elon,$edep);

}

# get_cmt_hdur_evdp
sub get_cmt_hdur_evdp {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my(undef,undef,$hdur)=split(" ",$cmt[3]);
  my(undef,$edep)=split(" ",$cmt[6]);

  return($hdur,$edep);
}


# get_cmt_moment : simply outputs M
sub get_cmt_moment {

  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my(undef,$Mrr)=split(" ",$cmt[7]);
  my(undef,$Mtt)=split(" ",$cmt[8]);
  my(undef,$Mpp)=split(" ",$cmt[9]);
  my(undef,$Mrt)=split(" ",$cmt[10]);
  my(undef,$Mrp)=split(" ",$cmt[11]);
  my(undef,$Mtp)=split(" ",$cmt[12]);

  return($Mrr,$Mtt,$Mpp,$Mrt,$Mrp,$Mtp);

}

# get_cmt_location: only outputs location
sub get_cmt_location {
 
  my ($cmt_file)=@_;

  open(CMT, "$cmt_file") or die("Error opening $cmt_file\n");
  @cmt = <CMT>;
  close(CMT);
  my(undef,$elat)=split(" ",$cmt[4]);
  my(undef,$elon)=split(" ",$cmt[5]);
  my(undef,$edep)=split(" ",$cmt[6]);

  return($elat,$elon,$edep);
}

# ($mw,$elat,$elon,$edep,$s1,$d1,$r1,$s2,$d2,$r2,$moment) = cmt2dc ($cmt_file)
sub cmt2dc{

  my ($cmt_file) = @_;
  my (@M,$elat,$elon,$edep,$mij2dc,$output);
  my ($s1,$d1,$r1,$s2,$d2,$r2,$moment,$mw);

  (@M[0..5],$elat,$elon,$edep) = get_cmt_short($cmt_file);
  $mij2dc = "${BIN}/mij2dc";
  if (not -f $mij2dc) {die("not mij2dc file\n");}
  $output = "mij2dc.tmp";
  open(MIJ,"|$mij2dc > $output");
  print MIJ "1 @M[0..5]\n";
  close(MIJ);
  (undef,$s1,$d1,$r1) = split(" ",`egrep output1 $output`);
  (undef,$s2,$d2,$r2) = split(" ",`egrep output2 $output`);
  (undef,undef,$moment) =
    split(" ",`egrep 'scalar moment' $output`);
  system("\\rm -f $output ");
  $mw = log($moment)/log(10)/1.5-10.73;
  return($mw,$elat,$elon,$edep,$s1,$d1,$r1,$s2,$d2,$r2,$moment);
}


# ($s1,$d1,$r1,$s2,$d2,$r2,$moment) = cmt2dc_simple($cmt_file)
sub cmt2dc_simple{

  my ($cmt_file) = @_;
  my (@M,$mij2dc,$output);
  my ($s1,$d1,$r1,$s2,$d2,$r2,$moment);

  (@M[0..5]) = get_cmt_moment($cmt_file);
  $mij2dc = "${BIN}/mij2dc";
  if (not -f $mij2dc) {die("no  $mij2dc file\n");}
  $output = "mij2dc.tmp";
  open(MIJ,"|$mij2dc > $output");
  print MIJ "1 @M[0..5]\n";
  close(MIJ);
  (undef,$s1,$d1,$r1) = split(" ",`egrep output1 $output`);
  (undef,$s2,$d2,$r2) = split(" ",`egrep output2 $output`);
  (undef,undef,$moment) =
    split(" ",`egrep 'scalar moment' $output`);
  system("\\rm -f $output ");
  return($s1,$d1,$r1,$s2,$d2,$r2,$moment);
}

# convert moment tensors to double couples
sub moment2dc{
  if (@_ != 6) {die("moment2dc(moment_tensor[0..5])\n");}
  my (@M) = @_;
  my ($s1,$d1,$r1,$s2,$d2,$r2,$moment,$mij2dc,$output);

  $mij2dc = "${BIN}/mij2dc";
  if (not -f $mij2dc) {die("no  $mij2dc file\n");}
  $output = "mij2dc.tmp";
  open(MIJ,"|$mij2dc > $output");
  print MIJ "1 @M[0..5]\n";
  close(MIJ);
  (undef,$s1,$d1,$r1) = split(" ",`egrep output1 $output`);
  (undef,$s2,$d2,$r2) = split(" ",`egrep output2 $output`);
  (undef,undef,$moment) =
    split(" ",`egrep 'scalar moment' $output`);
  system("\\rm -f $output ");
  return($s1,$d1,$r1,$s2,$d2,$r2,$moment);
}

# add tdiff time
sub tdiff{
  my ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tadd)=@_;
  $time = timegm($osec, $omin, $ohr, $oday , $omonth-1, $oyear);
  $time += ($tadd +$omsec/1000); #event_time in machine format
  $msec = sprintf("%03.0f",($time - (floor($time)))*1000);
  $time = floor($time);
  my($sec, $min, $hr, $day , $month, $year,$weekday,$jday) = gmtime($time);
  $month += 1;
  $year += 1900;
  $jday +=1;
  return ($year,$jday,$hr,$min,$sec,$msec);
}

# add tdiff time
sub tdiff2{
  my ($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tadd)=@_;
  $time = timegm($osec, $omin, $ohr, $oday , $omonth-1, $oyear);
  $time += ($tadd +$omsec/1000); #event_time in machine format
  $msec = sprintf("%03.0f",($time - (floor($time)))*1000);
  $time = floor($time);
  my($sec, $min, $hr, $day , $month, $year,$weekday,$jday) = gmtime($time);
  $month += 1;
  $year += 1900;
  $jday +=1;
  return ($year,$month,$day,$hr,$min,$sec,$msec);
}


sub mday2jday {
  my($oyear,$omonth,$oday)=@_;
  $omonth = $omonth-1;  #months range from 0..11
  $time_sec = timegm(3,3,3,$oday,$omonth,$oyear);
  @t = gmtime($time_sec);
  my($jday) = $t[7];
  $jday += 1;
  return ($jday);
}

sub jday2mday {
  my($year,$jday)=@_;
  $time = timegm(1,0,0,1,0,$year);
  $time=$time+($jday-1)*24*60*60;
  my @gmt_time=gmtime($time);
  return($gmt_time[4]+1,$gmt_time[3]);
}


1;
