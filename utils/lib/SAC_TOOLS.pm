package SAC_TOOLS;
# Qinya Liu, Caltech, May 2007

use warnings;
use Exporter;
use lib '/opt/seismo-util/lib/perl';
use SACLST;

@ISA = ('Exporter');
@EXPORT = ('write_station_psxy','write_station_pstext',
	   'match_data_and_syn','copyhdr_data_to_syn');


sub write_station_psxy{
  my(@sacfiles) = @_;
  my(@slat,@slon,$output);

  @slat = sac_files_get_values("stla",@sacfiles);
  @slon = sac_files_get_values("stlo",@sacfiles);

  $output = "";
  for ($i=0;$i<@slat;$i++) {
    $output .= "$slon[$i] $slat[$i]\n";
  }
  chomp($output);
  return($output);
}


sub write_station_pstext{
  my($options,$dxdy,@sacfiles)=@_;
  my(@slat,@slon,@stnm,$tlat,$tlon,$dx,$dy);

  @slat = sac_files_get_values("stla",@sacfiles);
  @slon = sac_files_get_values("stlo",@sacfiles);
  @stnm = sac_files_get_values("kstnm",@sacfiles);

  if (not $options) {$options = "10 0 4 CB";}

  ($dx,$dy) = split(" ",$dxdy);
  if (not defined $dx or not defined $dy ) {$dx = 0; $dy = 0;}

  $output = "";
  for ($i=0;$i<@slat;$i++) {
    $tlat = $slat[$i] + $dy; $tlon = $slon[$i] + $dx;
    $output .= "$tlon $tlat $options $stnm[$i]\n";
  }
  chomp($output); 
  return($output);
}

# @syn_files = match_data_and_syn("syn,","@data_files")

sub match_data_and_syn{
  my($syn_info,$data_files,$new_syn_info) = @_;
  my($syn_dir,$syn_suff,$data,$syn,$new_syn,@matched_data,@matched_syn,@matched_new_syn);
  @matched_data = ();
  @matched_syn = ();
  @matched_new_syn = ();

  ($syn_dir,$syn_suff) = split(/\,/,$syn_info);
  if (not $syn_dir) {$syn_dir = ".";}
  if ($syn_suff) {$syn_suff = ".$syn_suff";} else {$syn_suff = "";}
  if ($new_syn_info) {
       ($new_syn_dir,$new_syn_suff) = split(/\,/,$new_syn_info);
       if (not $new_syn_dir) {$new_syn_dir = ".";}
       if ($new_syn_suff) {$new_syn_suff = ".$new_syn_suff";} else {$new_syn_suff = "";}
   }

  (@data_files) = split(" ",$data_files);
  foreach $data (@data_files) {
    my($tmp) = sac_files_get_values("kstnm knetwk kcmpnm",$data);
    ($sta,$net,$comp) = split (" ",$tmp);
    $syn = "$syn_dir/$sta.$net.$comp${syn_suff}";
    if ($new_syn_info) {$new_syn = "$new_syn_dir/$sta.$net.$comp${new_syn_suff}";}
    if (-f $syn and -f $data) {
      if ($new_syn_info) { 
          if (-f $new_syn) {
            push(@matched_data,$data); push(@matched_syn,$syn); push(@matched_new_syn,$new_syn);}}
      else {
            push(@matched_data,$data); push(@matched_syn,$syn);}
    }
  }
  return("@matched_data","@matched_syn","@matched_new_syn");
}


sub copyhdr_data_to_syn{

  my($data_file,$syn_file,$headers) = @_;
  my(@headers,$header);

  (@headers) = split(" ",$headers);
  foreach $header (@headers) {
    if ($header !~ /^t\d+$/) {
      print "We could only change t? header at this point, not $header\n";
      return;}}
  open(SAC_SUB,"|sac2000 > /dev/null");
  print SAC_SUB "rh $data_file $syn_file\n";
  print SAC_SUB "copyhdr $headers from 1 \n wh\n";
  close(SAC_SUB);
}
